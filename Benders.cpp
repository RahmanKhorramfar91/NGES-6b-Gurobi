#include"Models_Funcs.h";

double Benders_Decomposition(GRBEnv* env)
{
	auto start = chrono::high_resolution_clock::now();
	vector<SP> Cuts;
	double MP_obj = 0; //MP objective in the current iterations
	int iter = 0;
	int iter_lim = 0; double gap = 1;
	double LB_obj = 0;
	double feas_obj = 0;

	MP_obj = Master_Problem(Cuts, env, 3600, 0.08, LB_obj, MP_obj);
	//bool org_UC_status = Setting::UC_active;
	//Setting::UC_active = false;
	//Integrated_Model(env, 0.06);
	//Setting::UC_active = org_UC_status;
	Setting::fix_some_E_NG_vars = true;
	feas_obj = Integrated_Model(env, 0.08);
	Setting::fix_some_E_NG_vars = false;
	EV::Benders_iter++;
	Setting::warm_start_active = true;

	//LB_obj = MP_init_heuristic(env);

	std::cout << endl;
	double gap0 = 0.06;

	int K = 1;
	double MP_gap = 0;
	while (true)
	{
		MP_obj = Master_Problem(Cuts, env, 3600, gap0, LB_obj, MP_gap);
		gap0 = std::max(0.01, gap0 * 0.8);
		if (iter_lim == K - 1)
		{
			gap0 = 0.01;
		}

		double Dual_stat = Dual_Subproblem(Cuts, env); // 0: unbounded, -1: feasible sol
		//double gap = std::abs(MP_obj0 - MP_obj1) / MP_obj1;

		
		if (iter_lim == K)
		{
			if (MP_gap > 0.0101)
			{
				return -1;
			}
			Setting::fix_some_E_NG_vars = true;
			feas_obj = Integrated_Model(env, Setting::cplex_gap);
			gap = (feas_obj - MP_obj) / feas_obj;
			Setting::fix_some_E_NG_vars = false;
			iter_lim = 0;
			K = std::max(1, int(K * 0.6));
		}
		EV::Benders_iter++;

		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		std::cout << " Elapsed time: " << Elapsed << "\t Iteration: " << iter << "\t MP Obj: " << MP_obj <<"\t UB: "<<feas_obj<< endl;
		if (gap < 0.01)
		{
			break;
		}
		if (Elapsed > Setting::CPU_limit)
		{
			EV::MIP_gap = gap;
			break;
		}

		iter++; iter_lim++;
	}

	return feas_obj;
}


double Master_Problem(vector<SP> Cuts, GRBEnv* env, double CPU_time, double gap0, double LB_obj, double& MP_gap)
{
	//auto start = chrono::high_resolution_clock::now();
#pragma region Fetch Data
	std::map<string, int> sym2pltType = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7},
		{"CT",8},{"CC",9},{"CC-CCS",10},{"solar-UPV",11},{"wind-new",12},
		{"wind-offshore-new",13},{"hydro-new",14},{"nuclear-new",15} };
	std::map<int, string> pltType2sym = { {0,"ng"},{1,"dfo"},
{2,"solar"},{3,"wind"},{4,"wind_offshore"},{5,"hydro"},{6,"coal"},{7,"nuclear"} };
	vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Params::Enodes.size();
	int nPlt = (int)Params::Plants.size();
	int nBr = (int)Params::Branches.size();
	int neSt = (int)Params::Estorage.size();
	vector<int> Tg = Params::Tg;
	vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	double pi = 3.141592;
	int nGnode = (int)Params::Gnodes.size();
	double WACC = Params::WACC;
	int trans_unit_cost = Params::trans_unit_cost;
	int trans_line_lifespan = Params::trans_line_lifespan;
	double NG_price = Params::NG_price;
	double dfo_pric = Params::dfo_pric;
	double coal_price = Params::coal_price;
	double nuclear_price = Params::nuclear_price;
	double E_curt_cost = Params::E_curt_cost;
	double G_curt_cost = Params::G_curt_cost;
	double pipe_per_mile = Params::pipe_per_mile;
	int pipe_lifespan = Params::pipe_lifespan;
	int battery_lifetime = Params::battery_lifetime;
	map<int, vector<int>> Le = Params::Le;
	vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion


	//GRBEnv* env = 0;

	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_NGobj(0);
	GRBLinExpr exp_Eobj(0);

	Populate_EV_SP(Model);
	Populate_GV(Model);
	NG_Module(Model, exp_NGobj);

	// coupling constraints
	GRBLinExpr ex_xi(0);
	GRBLinExpr ex_NG_emis(0);
	GRBLinExpr ex_E_emis(0);
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);

	//if (Setting::Case == 1)
	//{
	//	Model.addConstr(CV::xi == 100); // just to populate CV::xi
	//} else if to next 


	Model.addConstr(CV::xi == ex_xi);



	if (Setting::Case == 3)
	{// the original model
		Model.addConstr(ex_E_emis + ex_NG_emis <=(1-Setting::Emis_redu_goal) * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}

#pragma region Set some variables
	// 1) existing types can not be established because there are new equivalent types
	// 2) new types cannot be decommissioned
	// 3) no new wind-offshore on some location (not yet implemented)
	// 4) no new nuclear
	// 5) turn off df0/coal/wind-offshore. The data is not correct on wind-offshore, there is no existing wind-offshore in the region
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			//// do not allow establishment off-shore wind and hydro plants
			//Model.addConstr(EV::Xdec[n][5] == 0);
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore")
			{
				Model.addConstr(EV::Xop[n][i] == 0);
			}

			if (Plants[i].type == "hydro-new")
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}

			if (Plants[i].type == "wind-offshore-new")
			{
				if (n == 0 || n == 2)// except for Massachusetts and Connecticus
				{
					continue;
				}
				Model.addConstr(EV::Xest[n][i] == 0);
			}

			/*if (Plants[i].type == "nuclear-new")
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}*/

			if (Plants[i].is_exis == 1)
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}
			else
			{
				Model.addConstr(EV::Xdec[n][i] == 0);
			}
		}
	}

	//if (Setting::MP_init_heuristic)
	//{
	//	for (int n = 0; n < nEnode; n++)
	//	{
	//		for (int i = 0; i < nPlt; i++)
	//		{
	//			Model.addConstr(EV::Xest[n][i] == EV::val_Xest[n][i]);
	//			Model.addConstr(EV::Xdec[n][i] == EV::val_Xdec[n][i]);
	//			Model.addConstr(EV::Xop[n][i] == EV::val_Xop[n][i]);
	//			for (int t = 0; t < Te.size(); t++)
	//			{
	//				Model.addConstr(EV::X[n][t][i] == EV::val_X[n][t][i]);
	//				Model.addConstr(EV::Xup[n][t][i] == EV::val_Xup[n][t][i]);
	//				Model.addConstr(EV::Xdown[n][t][i] == EV::val_Xdown[n][t][i]);
	//			}
	//		}
	//	}
	//}
	//Setting::MP_init_heuristic = false;//to only apply in the first iteration
#pragma endregion


#pragma region Objective Function
	GRBLinExpr ex_est(0);
	GRBLinExpr ex_decom(0);
	GRBLinExpr ex_fix(0);
	GRBLinExpr ex_startup(0);
	GRBLinExpr ex_var(0);
	GRBLinExpr ex_thermal_fuel(0);
	GRBLinExpr ex_emis(0);
	GRBLinExpr ex_shedd(0);
	GRBLinExpr ex_trans(0);
	GRBLinExpr ex_elec_str(0);


	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			// Investment/decommission cost of plants
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore") {
				continue;
			}
			double s1 = (double)std::pow(1.0 / (1 + WACC), Plants[i].lifetime);
			double capco = WACC / (1 - s1) * Plants[i].Reg_coeffs_per_state[n];
			ex_est += capco * Plants[i].capex * EV::Xest[n][i];
			ex_decom += Plants[i].decom_cost * EV::Xdec[n][i];
		}

		// fixed cost (annual, so no iteration over time)
		for (int i = 0; i < nPlt; i++)
		{
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore") {
				continue;
			}
			ex_fix += Plants[i].fix_cost * EV::Xop[n][i];
		}
		// var+fuel costs of plants
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
					Plants[i].type == "wind_offshore") {
					continue;
				}
				// var cost
				ex_var += time_weight[t] * Plants[i].var_cost * EV::prod[n][t][i];

				// fuel price to be updated later (dollar per thousand cubic feet=MMBTu)
				// NG fuel is handle in the NG network for case 2 and 3
				//if (Plants[i].type == "dfo")
				//{
				//	ex_thermal_fuel += time_weight[t] * dfo_pric * Plants[i].heat_rate * EV::prod[n][t][i];

				//}
				//else if (Plants[i].type == "coal")
				//{
				//	ex_thermal_fuel += time_weight[t] * coal_price * Plants[i].heat_rate * EV::prod[n][t][i];
				//}
				if (Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
				{
					ex_thermal_fuel += time_weight[t] * nuclear_price * Plants[i].heat_rate * EV::prod[n][t][i];
				}
				else if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
				{
					if (Setting::Case == 1)
					{
						ex_thermal_fuel += time_weight[t] * NG_price * Plants[i].heat_rate * EV::prod[n][t][i];
					}
				}

				if (Plants[i].type == "ng" || Plants[i].type == "CT" ||
					Plants[i].type == "CC" || Plants[i].type == "CC-CCS" ||
					Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
				{
					// startup cost
					ex_startup += time_weight[t] * Plants[i].startup_cost * EV::Xup[n][t][i];
				}

			}

			// load curtailment cost
			ex_shedd += time_weight[t] * E_curt_cost * EV::curtE[n][t];
		}


		// storage cost
		for (int r = 0; r < neSt; r++)
		{
			double s1 = (double)std::pow(1.0 / (1 + WACC), battery_lifetime);
			double capco = WACC / (1 - s1);
			ex_elec_str += capco * Estorage[r].power_cost * EV::YeCD[n][r] + capco * Estorage[r].energy_cost * EV::YeLev[n][r];
			ex_elec_str += Estorage[r].pFOM * EV::YeCD[n][r] + Estorage[r].eFOM * EV::YeLev[n][r]; // fixed cost per kw per year
		}
	}

	// investment cost of transmission lines
	for (int b = 0; b < Branches.size(); b++)
	{
		double s1 = (double)std::pow(1.0 / (1 + WACC), trans_line_lifespan);
		double capco = WACC / (1 - s1);
		int fb = Branches[b].from_bus;
		int tb = Branches[b].to_bus;
		// find the index of tb in the adj_buses of the fb
		int tbi = std::find(Enodes[fb].adj_buses.begin(), Enodes[fb].adj_buses.end(), tb) - Enodes[fb].adj_buses.begin();
		ex_trans += capco * trans_unit_cost * Branches[b].maxFlow * Branches[b].length * EV::Ze[b];
	}
	exp_Eobj = ex_est + ex_decom + ex_fix + ex_startup + ex_emis + ex_var + ex_thermal_fuel + ex_shedd + ex_trans + ex_elec_str;

	Model.addConstr(EV::est_cost == ex_est);
	Model.addConstr(EV::est_trans_cost == ex_trans);
	Model.addConstr(EV::decom_cost == ex_decom);
	Model.addConstr(EV::fixed_cost == ex_fix);
	Model.addConstr(EV::var_cost == ex_var);
	Model.addConstr(EV::thermal_fuel_cost == ex_thermal_fuel);
	Model.addConstr(EV::dfo_coal_emis_cost == ex_emis);
	Model.addConstr(EV::shedding_cost == ex_shedd);
	Model.addConstr(EV::elec_storage_cost == ex_elec_str);
	Model.addConstr(EV::startup_cost == ex_startup);

	Model.addConstr(exp_Eobj == EV::e_system_cost);
	Model.setObjective(exp_Eobj + exp_NGobj, GRB_MINIMIZE);

#pragma endregion

#pragma region MP Electricity local constraints
	// C1: number of generation units at each node	
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			// is plant type i part of initial plants types of node n:
			bool isInd = std::find(Enodes[n].Init_plt_ind.begin(), Enodes[n].Init_plt_ind.end(), i) != Enodes[n].Init_plt_ind.end();
			if (isInd)
			{
				// find the index of plant in the set of Init_plants of the node
				string tp = pltType2sym[i];
				int ind1 = std::find(Enodes[n].Init_plt_types.begin(), Enodes[n].Init_plt_types.end(), tp) - Enodes[n].Init_plt_types.begin();
				Model.addConstr(EV::Xop[n][i] == Enodes[n].Init_plt_count[ind1] - EV::Xdec[n][i] + EV::Xest[n][i]);
			}
			else
			{
				Model.addConstr(EV::Xop[n][i] == -EV::Xdec[n][i] + EV::Xest[n][i]);
			}

			for (int t = 0; t < Te.size(); t++)
			{
				// (w/o UC version) keep this to avoid power generation from hydro-new, dfor, etc.
				Model.addConstr(EV::prod[n][t][i] <= Plants[i].Pmax * EV::Xop[n][i]);
				if (!Setting::UC_active)
				{
					//(w/o UC version)
					if (t > 0)
					{
						//Model.addConstr(Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i] >= -EV::prod[n][t][i] + EV::prod[n][t - 1][i]);
						Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
						Model.addConstr(-EV::prod[n][t][i] + EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
					}
				}
			}
		}
	}

	if (Setting::UC_active)
	{
		//C2, C3, C4, C5: UC,  production limit, ramping for thermal units (ng, CT, CC, CC-CCS, nuclear)	
		for (int n = 0; n < nEnode; n++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				for (int t = 0; t < Te.size(); t++)
				{
					if (Plants[i].type == "ng" || Plants[i].type == "CT" ||
						Plants[i].type == "CC" || Plants[i].type == "CC-CCS" ||
						Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
					{
						// UC 1
						if (t > 0)
						{
							Model.addConstr(EV::X[n][t][i] - EV::X[n][t - 1][i] == EV::Xup[n][t][i] - EV::Xdown[n][t][i]);
						}
						// UC 2
						Model.addConstr(EV::X[n][t][i] <= EV::Xop[n][i]);
						// output limit
						Model.addConstr(EV::prod[n][t][i] >= Plants[i].Pmin * Plants[i].Pmax * EV::X[n][t][i]);
						Model.addConstr(EV::prod[n][t][i] <= Plants[i].Pmax * EV::X[n][t][i]);
						//ramping
						if (t > 0)
						{
							double mx1 = std::max(Plants[i].Pmin, Plants[i].rampU);
							Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * (EV::X[n][t][i] - EV::Xup[n][t][i]) + mx1 * Plants[i].Pmax * EV::Xup[n][t][i]);
							Model.addConstr(-EV::prod[n][t][i] + EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * (EV::X[n][t][i] - EV::Xup[n][t][i]) + mx1 * Plants[i].Pmax * EV::Xup[n][t][i]);
						}
					}
				}
			}
		}
	}


	////C3, C4: production limit, ramping	
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		for (int i = 0; i < nPlt; i++)
	//		{
	//			//Model.addConstr(EV::prod[n][t][i] >= Plants[i].Pmin * EV::Xop[n][i]); since we don't consider unit commitment in this model
	//			Model.addConstr(EV::prod[n][t][i] <= Plants[i].Pmax * EV::Xop[n][i]);

	//			//if (t > 0 && !Setting::heuristics1_active) 
	//			if (t > 0)
	//			{
	//				//Model.addConstr(Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i] >= -EV::prod[n][t][i] + EV::prod[n][t - 1][i]);
	//				Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
	//				Model.addConstr(-EV::prod[n][t][i] + EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
	//			}
	//		}
	//	}
	//}

	//// C 3.5: max yearly generation for existing plants
	//for (int i = 0; i < nPlt; i++)
	//{
	//	GRBLinExpr ex0(0);
	//	if (Plants[i].is_exis != 1) { continue; }

	//	for (int n = 0; n < nEnode; n++)
	//	{
	//		for (int t = 0; t < Te.size(); t++)
	//		{
	//			ex0 += time_weight[t] * EV::prod[n][t][i];
	//		}
	//	}
	//	Model.addConstr(ex0 <= Plants[i].max_yearly_gen);
	//}

	// C5, C6: flow limit for electricity
	//C7: power balance
	// C8: flow equation
	// C9: phase angle (theta) limits. already applied in the definition of the variable	
	// C10: VRE production profile
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "solar" || Plants[i].type == "solar-UPV")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::Xop[n][i]);
				}
				if (Plants[i].type == "wind" || Plants[i].type == "wind-new")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::Xop[n][i]);
				}
				if (Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][0] * Plants[i].Pmax * EV::Xop[n][i]);
				}
				/*if (Plants[i].type == "hydro"|| Plants[i].type == "hydro-new")
				{
					Model.addConstr(EV::prod[n][t][i] <=Plants[i].Pmax * EV::Xop[n][i]);
				}*/
			}
		}
	}


	// C11: demand curtainlment constraint
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			double dem = Enodes[n].demand[Te[t]];
			//dem = dem/1.896;
			Model.addConstr(EV::curtE[n][t] <= dem);
		}
	}

	//C12: RPS constraints
	GRBLinExpr exp3(0);
	double total_demand = 0;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			total_demand += time_weight[t] * Enodes[n].demand[Te[t]];

			for (int i = 0; i < nPlt; i++)
			{
				//exp2 += time_weight[t] * Plants[i].emis_rate * EV::prod[n][t][i];
				if (Plants[i].type == "solar" || Plants[i].type == "wind" ||
					Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new" ||
					Plants[i].type == "hydro" || Plants[i].type == "solar-UPV" ||
					Plants[i].type == "wind-new" || Plants[i].type == "hydro-new")
				{
					exp3 += time_weight[t] * EV::prod[n][t][i]; // production from VRE
				}
			}
		}
	}
	Model.addConstr(exp3 >= Setting::RPS * total_demand);

	// C14,C15,C16 storage constraints
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int r = 0; r < neSt; r++)
			{
				/*if (t == Te.size() - 1)
				{
					Model.addConstr(EV::eSlev[n][t][r] <= EV::eSlev[n][0][r]);
				}*/
				if (t == 0)
				{
					Model.addConstr(EV::eSlev[n][t][r] ==
						Estorage[r].eff_ch * EV::eSch[n][t][r] -
						(EV::eSdis[n][t][r] / Estorage[r].eff_disCh));
				}
				else
				{
					Model.addConstr(EV::eSlev[n][t][r] == EV::eSlev[n][t - 1][r] +
						Estorage[r].eff_ch * EV::eSch[n][t][r] -
						(EV::eSdis[n][t][r] / Estorage[r].eff_disCh));
				}
				Model.addConstr(EV::YeCD[n][r] >= EV::eSdis[n][t][r]);
				Model.addConstr(EV::YeCD[n][r] >= EV::eSch[n][t][r]);
				Model.addConstr(EV::YeLev[n][r] >= EV::eSlev[n][t][r]);
			}
		}
	}

	// start and ending of storage should be the same in case of using rep. days
	for (int n = 0; n < nEnode; n++)
	{
		for (int r = 0; r < neSt; r++)
		{
			for (int k = 0; k < Tg.size(); k++)
			{
				int str = k * 24;
				int end = (k + 1) * 24 - 1;
				Model.addConstr(EV::eSlev[n][str][r] == EV::eSlev[n][end][r]);
			}
		}
	}


	// C17: if an storage established or not
	for (int n = 0; n < nEnode; n++)
	{
		for (int r = 0; r < neSt; r++)
		{
			Model.addConstr(EV::YeLev[n][r] <= 10e8 * EV::YeStr[n][r]);
			Model.addConstr(EV::YeCD[n][r] <= EV::YeLev[n][r]);
		}
	}
#pragma endregion

#pragma region VIs
	for (int t = 0; t < Te.size(); t++)
	{
		GRBLinExpr ex_vi(0);
		double dem = 0;
		for (int n = 0; n < nEnode; n++)
		{
			//GRBLinExpr exp_prod(0);
			for (int i = 0; i < nPlt; i++) { ex_vi += EV::prod[n][t][i]; }

			//GRBLinExpr ex_store(0);
			for (int r = 0; r < neSt; r++) { ex_vi += EV::eSdis[n][t][r] - EV::eSch[n][t][r]; }

			dem += Enodes[n].demand[Te[t]];
			ex_vi += EV::curtE[n][t];
		}
		Model.addConstr(ex_vi == dem);
	}

	Model.addConstr(exp_Eobj + exp_NGobj >= LB_obj);
#pragma endregion

#pragma region Add cuts from SPs
	// define the contribution from subproblem
	//GRBVar Chi = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);

	//Chi = Model.addVars(2, GRB_CONTINUOUS);
	//Model.addConstr(Chi[0] >= ex_cut);
	//Chi[0] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);

	//bool Multi_cut = false;  // multicut implementation or single cut
	GRBLinExpr ex_cut(0);
	for (int c = 0; c < Cuts.size(); c++)
	{
		/*if (c==1)
		{
			continue;
		}*/
		for (int t = 0; t < Te.size(); t++)
		{
			for (int n = 0; n < nEnode; n++)
			{
				GRBLinExpr exp_prod(0);
				for (int i = 0; i < nPlt; i++)
				{
					exp_prod += EV::prod[n][t][i];
				}
				double dem = Enodes[n].demand[Te[t]];
				GRBLinExpr str(0);
				for (int r = 0; r < neSt; r++)
				{
					str += EV::eSdis[n][t][r] - EV::eSch[n][t][r];
				}
				if (Cuts[c].dual_val_theta1[n][t] > 0)
				{
					ex_cut += (dem - str - exp_prod - EV::curtE[n][t]) * Cuts[c].dual_val_theta1[n][t];
				}

				//lhs += dem;
				if (Cuts[c].dual_val_theta2[n][t] > 0)
				{
					ex_cut -= (dem - str - exp_prod - EV::curtE[n][t]) * Cuts[c].dual_val_theta2[n][t];
				}
				//lhs += dem;
				ex_cut -= pi * Cuts[c].dual_val_eta1[n][t];
				ex_cut -= pi * Cuts[c].dual_val_eta2[n][t];
			}

			for (int b = 0; b < nBr; b++)
			{

				if (Branches[b].is_exist == 1)
				{
					ex_cut -= Branches[b].maxFlow * Cuts[c].dual_val_delta11[b][t];
					ex_cut -= Branches[b].maxFlow * Cuts[c].dual_val_delta12[b][t];
					//lhs += Branches[b].maxFlow * Cuts[c].dual_val_delta11[b][t];
					//lhs += Branches[b].maxFlow * Cuts[c].dual_val_delta12[b][t];
				}
				else
				{
					double Big_M = std::abs(Branches[b].suscep * 24);
					ex_cut -= Branches[b].maxFlow * EV::Ze[b] * Cuts[c].dual_val_delta21[b][t];
					ex_cut -= Branches[b].maxFlow * EV::Ze[b] * Cuts[c].dual_val_delta22[b][t];
					ex_cut -= Big_M * (1 - EV::Ze[b]) * Cuts[c].dual_val_zeta11[b][t];
					ex_cut -= Big_M * (1 - EV::Ze[b]) * Cuts[c].dual_val_zeta12[b][t];
				}
			}
			if (Setting::multi_cut_active == true)
			{
				Model.addConstr(ex_cut <= 0); // all cuts are feasibility cuts
				ex_cut.clear();
			}
		}
		if (Setting::multi_cut_active == false)
		{
			Model.addConstr(ex_cut <= 0); // all cuts are feasibility cuts
			ex_cut.clear();
		}
	}

#pragma endregion


#pragma region Add warm start solution
	if (Setting::warm_start_active)
	{
		for (int n = 0; n < nEnode; n++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				//Model.set(GRB_DoubleAttr_Start, EV::Xop[n],EV::val_Xop[n],nPlt);
				EV::Xop[n][i].set(GRB_DoubleAttr_Start, EV::val_Xop[n][i]);
				EV::Xest[n][i].set(GRB_DoubleAttr_Start, EV::val_Xest[n][i]);
				EV::Xdec[n][i].set(GRB_DoubleAttr_Start, EV::val_Xdec[n][i]);
				for (int t = 0; t < Te.size(); t++)
				{
					//EV::prod[n][t][i].set(GRB_DoubleAttr_Start, EV::val_prod[n][t][i]);
					EV::X[n][t][i].set(GRB_DoubleAttr_Start, EV::val_X[n][t][i]);
				}
			}
		}
		for (int b = 0; b < nBr; b++)
		{
			EV::Ze[b].set(GRB_DoubleAttr_Start, EV::val_Ze[b]);
			GV::Zg[b].set(GRB_DoubleAttr_Start, GV::val_Zg[b]);
		}
		Model.update();
	}
#pragma endregion


#pragma region Solve MP
	Model.set(GRB_IntParam_OutputFlag, 0);
	Model.set(GRB_DoubleParam_TimeLimit, CPU_time);
	Model.set(GRB_DoubleParam_MIPGap, gap0);
	Model.optimize();
	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		//std::cout << "Failed to optimize IM!!!" << endl;
		//std::cout << Model.get(GRB_IntAttr_Status);
	}

	double obj_val = Model.get(GRB_DoubleAttr_ObjVal);

	if (!Setting::relax_int_vars)
	{
		MP_gap = Model.get(GRB_DoubleAttr_MIPGap);
	}

	if (MP_gap > 0.1)
	{
		return -1;
	}
	//auto end = chrono::high_resolution_clock::now();
	//double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	//std::cout << "Elapsed time: " << Elapsed << endl;
	//std::cout << "\t MP Obj Value:" << obj_val << endl;
	//std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;
	////std::cout << "\t Chi value = " << Chi.get(GRB_DoubleAttr_X) << endl;
	//std::cout << "\t NG emission = " << CV::NG_emis.get(GRB_DoubleAttr_X) << endl;
	//std::cout << "\t E emission = " << CV::E_emis.get(GRB_DoubleAttr_X) << endl;
	//std::cout << "\t Trans. Est. Cost = " << EV::est_trans_cost.get(GRB_DoubleAttr_X) << endl;
#pragma endregion

	//Get_EV_vals(Model);
	//Get_GV_vals(Model);
	/*Print_Results(0, obj_val);*/

#pragma region Get the values of MP variables	
	EV::val_Xop = new double* [nEnode];
	EV::val_Xest = new double* [nEnode];
	EV::val_Xdec = new double* [nEnode];
	EV::val_X = new double** [nEnode];
	EV::val_Xup = new double** [nEnode];
	EV::val_Xdown = new double** [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_Xop[n] = new double[nPlt]();
		EV::val_Xest[n] = new double[nPlt]();
		EV::val_Xdec[n] = new double[nPlt]();
		for (int i = 0; i < nPlt; i++)
		{
			EV::val_Xop[n][i] = std::round(EV::Xop[n][i].get(GRB_DoubleAttr_X));
			EV::val_Xest[n][i] = std::round(EV::Xest[n][i].get(GRB_DoubleAttr_X));
			EV::val_Xdec[n][i] = std::round(EV::Xdec[n][i].get(GRB_DoubleAttr_X));
		}
	}
	EV::val_Ze = new double[nBr]();
	//EV::val_num_est_trans = 0;
	for (int b = 0; b < nBr; b++)
	{
		EV::val_Ze[b] = std::round(EV::Ze[b].get(GRB_DoubleAttr_X));
	}

	EV::val_prod = new double** [nEnode];
	double total_yearly_prod = 0;
	EV::val_total_prod = new double[nPlt]();
	double flowGE = 0;
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_prod[n] = new double* [Te.size()];
		EV::val_X[n] = new double* [Te.size()];
		EV::val_Xup[n] = new double* [Te.size()];
		EV::val_Xdown[n] = new double* [Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			//if (t > periods2print) { break; }

			EV::val_prod[n][t] = new double[nPlt]();
			EV::val_X[n][t] = new double[nPlt]();
			EV::val_Xup[n][t] = new double[nPlt]();
			EV::val_Xdown[n][t] = new double[nPlt]();
			for (int i = 0; i < nPlt; i++)
			{
				EV::val_total_prod[i] += time_weight[t] * EV::prod[n][t][i].get(GRB_DoubleAttr_X);
				total_yearly_prod += time_weight[t] * EV::prod[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_prod[n][t][i] = EV::prod[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_X[n][t][i] = EV::X[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_Xup[n][t][i] = EV::Xup[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_Xdown[n][t][i] = EV::Xdown[n][t][i].get(GRB_DoubleAttr_X);
			}
		}
	}
	/*for (int i = 0; i < nPlt; i++)
	{
		if (EV::val_total_prod[i] > 10)
		{
			cout << Plants[i].type << ": \t " << EV::val_total_prod[i] << endl;
		}
	}*/
	//cout << "\t\t total yearly product: " << total_yearly_prod << endl;

	EV::val_curtE = new double* [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_curtE[n] = new double[Te.size()]();
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_total_curt += time_weight[t] * EV::curtE[n][t].get(GRB_DoubleAttr_X);
			EV::val_curtE[n][t] = EV::curtE[n][t].get(GRB_DoubleAttr_X);
		}
	}

	EV::val_eSch = new double** [nEnode];
	EV::val_eSdis = new double** [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_eSch[n] = new double* [Te.size()];
		EV::val_eSdis[n] = new double* [Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_eSch[n][t] = new double[neSt]();
			EV::val_eSdis[n][t] = new double[neSt]();
			for (int i = 0; i < neSt; i++)
			{
				EV::val_eSch[n][t][i] = EV::eSch[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_eSdis[n][t][i] = EV::eSdis[n][t][i].get(GRB_DoubleAttr_X);
			}
		}
	}



	GV::val_supply = new double* [nGnode];
	for (int k = 0; k < nGnode; k++)
	{
		GV::val_supply[k] = new double[Tg.size()]();
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			GV::val_supply[k][tau] = GV::supply[k][tau].get(GRB_DoubleAttr_X);
		}
	}

	GV::val_Zg = new double[Params::PipeLines.size()]();
	for (int i = 0; i < Params::PipeLines.size(); i++)
	{
		GV::val_num_est_pipe += GV::Zg[i].get(GRB_DoubleAttr_X);
	}

	/*
		GV::val_flowGE = new double** [nGnode];
		for (int k = 0; k < nGnode; k++)
		{
			GV::val_flowGE[k] = new double* [nEnode];
			for (int kp : Gnodes[k].adjE)
			{
				GV::val_flowGE[k][kp] = new double[Tg.size()]();
				for (int tau = 0; tau < Tg.size(); tau++)
				{
					GV::val_flowGE[k][kp][tau] = GV::flowGE[k][kp][tau].get(GRB_DoubleAttr_X);
				}
			}
		}*/
#pragma endregion


#pragma region Check Chi val
		//double Chi_val = 0;

		//for (int c = 0; c < Cuts.size(); c++)
		//{
		//	//for (int i = 0; i < nPlt; i++)
		//	//{
		//	//	if (Plants[i].is_exis != 1) { continue; } // only applies to existing plants
		//	//	Chi_val -= Plants[i].max_yearly_gen * Cuts[c].dual_val_alpha[i];
		//	//}

		//	for (int n = 0; n < nEnode; n++)
		//	{
		//		for (int t = 0; t < Te.size(); t++)
		//		{
		//			double exp_prod = 0;
		//			for (int i = 0; i < nPlt; i++)
		//			{
		//				exp_prod += EV::val_prod[n][t][i];
		//			}
		//			double dem = Enodes[n].demand[Te[t]];
		//			double str = 0;
		//			for (int r = 0; r < neSt; r++)
		//			{
		//				str += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r];
		//			}
		//			Chi_val += (dem - str - exp_prod - EV::val_curtE[n][t]) * Cuts[c].dual_val_theta1[n][t];
		//			Chi_val -= (dem - str - exp_prod - EV::val_curtE[n][t]) * Cuts[c].dual_val_theta2[n][t];
		//			Chi_val -= pi * Cuts[c].dual_val_eta1[n][t];
		//			Chi_val -= pi * Cuts[c].dual_val_eta2[n][t];
		//		}
		//	}

		//	for (int b = 0; b < nBr; b++)
		//	{
		//		for (int t = 0; t < Te.size(); t++)
		//		{
		//			if (Branches[b].is_exist == 1)
		//			{
		//				Chi_val -= Branches[b].maxFlow * Cuts[c].dual_val_delta11[b][t];
		//				Chi_val -= Branches[b].maxFlow * Cuts[c].dual_val_delta12[b][t];
		//			}
		//			else
		//			{
		//				double Big_M = std::abs(Branches[b].suscep * 24);
		//				Chi_val -= Branches[b].maxFlow * EV::val_Ze[b] * Cuts[c].dual_val_delta21[b][t];
		//				Chi_val -= Branches[b].maxFlow * EV::val_Ze[b] * Cuts[c].dual_val_delta22[b][t];
		//				Chi_val += Big_M * (1 - EV::val_Ze[b]) * Cuts[c].dual_val_zeta11[b][t];
		//				Chi_val -= Big_M * (1 - EV::val_Ze[b]) * Cuts[c].dual_val_zeta12[b][t];
		//			}
		//		}
		//	}
		//}
		//std::cout << "Chi rhs: " << Chi_val << endl;
#pragma endregion

	return obj_val;
}


double Dual_Subproblem(vector<SP>& Cuts, GRBEnv* env)
{

	SP sp_new;
	double eps = 0.02;
	auto start = chrono::high_resolution_clock::now();

#pragma region Fetch Data
	vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Enodes.size();
	int nPlt = (int)Plants.size();
	int nBr = (int)Branches.size();
	int neSt = (int)Estorage.size();
	vector<int> Tg = Params::Tg;
	vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	const double pi = 3.1415926535897;
	int nGnode = (int)Gnodes.size();
	double WACC = Params::WACC;
	int trans_unit_cost = Params::trans_unit_cost;
	int trans_line_lifespan = Params::trans_line_lifespan;
	double NG_price = Params::NG_price;
	double dfo_pric = Params::dfo_pric;
	double coal_price = Params::coal_price;
	double E_curt_cost = Params::E_curt_cost;
	double G_curt_cost = Params::G_curt_cost;
	double pipe_per_mile = Params::pipe_per_mile;
	int pipe_lifespan = Params::pipe_lifespan;
	map<int, vector<int>> Le = Params::Le;
	vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion
	/*GRBEnv* env = 0;
	env = new GRBEnv();*/
	GRBModel Model = GRBModel(env);
#pragma region Define Decision Variables
	/*SP::alpha = new GRBVar[nPlt];
	SP::beta = new GRBVar * *[nEnode];
	SP::gamma1 = new GRBVar * *[nEnode];
	SP::gamma2 = new GRBVar * *[nEnode];*/
	SP::delta11 = new GRBVar * [nBr];
	SP::delta12 = new GRBVar * [nBr];
	SP::delta21 = new GRBVar * [nBr];
	SP::delta22 = new GRBVar * [nBr];
	SP::theta1 = new GRBVar * [nEnode];
	SP::theta2 = new GRBVar * [nEnode];
	SP::zeta11 = new GRBVar * [nBr];
	SP::zeta21 = new GRBVar * [nBr];
	SP::zeta12 = new GRBVar * [nBr];
	SP::zeta22 = new GRBVar * [nBr];
	SP::eta1 = new GRBVar * [nEnode];
	SP::eta2 = new GRBVar * [nEnode];
	SP::eta3 = new GRBVar[Te.size()];
	/*SP::pi = new GRBVar * *[nEnode];
	SP::phi = new GRBVar * [nEnode];
	SP::rho = new GRBVar * *[nGnode];
	SP::omega = Model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	SP::tau = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);*/

	// free variable
	for (int t = 0; t < Te.size(); t++)
	{
		SP::eta3[t] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
	}

	for (int n = 0; n < nEnode; n++)
	{
		/*SP::beta[n] = new GRBVar * [Te.size()];
		SP::gamma1[n] = new GRBVar * [Te.size()]; SP::gamma2[n] = new GRBVar * [Te.size()];*/
		SP::theta1[n] = new GRBVar[Te.size()];
		SP::theta2[n] = new GRBVar[Te.size()];
		SP::eta1[n] = new GRBVar[Te.size()];
		SP::eta2[n] = new GRBVar[Te.size()];



		//SP::pi[n] = new GRBVar * [Te.size()];
		//SP::phi[n] = new GRBVar[Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			//SP::beta[n][t] = new GRBVar[nPlt];
			//SP::gamma1[n][t] = new GRBVar[nPlt]; SP::gamma2[n][t] = new GRBVar[nPlt];

			// free dual variable
			SP::theta1[n][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);

			SP::theta2[n][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);

			SP::eta1[n][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::eta2[n][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			////SP::pi[n][t] = new GRBVar[nPlt];
			//SP::phi[n][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			//for (int i = 0; i < nPlt; i++)
			//{
			//	SP::beta[n][t][i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			//	SP::gamma1[n][t][i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			//	SP::gamma2[n][t][i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			//	SP::pi[n][t][i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			//}
		}
	}


	for (int b = 0; b < nBr; b++)
	{
		SP::delta11[b] = new GRBVar[Te.size()];
		SP::delta12[b] = new GRBVar[Te.size()];
		SP::delta21[b] = new GRBVar[Te.size()];
		SP::delta22[b] = new GRBVar[Te.size()];
		SP::zeta11[b] = new GRBVar[Te.size()];
		SP::zeta12[b] = new GRBVar[Te.size()];
		SP::zeta21[b] = new GRBVar[Te.size()];
		SP::zeta22[b] = new GRBVar[Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			SP::delta11[b][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::delta12[b][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::delta21[b][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::delta22[b][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);

			SP::zeta11[b][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::zeta12[b][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::zeta21[b][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::zeta22[b][t] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		}
	}

	//for (int k = 0; k < nGnode; k++)
	//{
	//	SP::rho[k] = new GRBVar * [nEnode];
	//	for (int n = 0; n < nEnode; n++)
	//	{
	//		SP::rho[k][n] = new GRBVar[Tg.size()];
	//		for (int tau = 0; tau < Tg.size(); tau++)
	//		{
	//			// free variable
	//			SP::rho[k][n][tau] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
	//		}
	//	}
	//}
#pragma endregion



#pragma region Set some variable (used in verification step)
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			//Model.addConstr(SP::delta21[b][t] == 0);
			//Model.addConstr(SP::delta22[b][t] == 0);
			//Model.addConstr(SP::delta11[b][t] == 0);
			//Model.addConstr(SP::delta12[b][t] == 0);

			//Model.addConstr(SP::zeta11[b][t] == 0);
			//Model.addConstr(SP::zeta12[b][t] == 0);
			//Model.addConstr(SP::zeta21[b][t] == 0);
			//Model.addConstr(SP::zeta22[b][t] == 0);
		}
	}

	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			//Model.addConstr(SP::theta[n][t] == SP::dual_val_theta[n][t]);
			//Model.addConstr(SP::phi[n][t] == SP::dual_val_phi[n][t]);
			for (int i = 0; i < nPlt; i++)
			{
				//Model.addConstr(SP::pi[n][t][i] == 0);
				//Model.addConstr(SP::gamma1[n][t][i] == 0);
				//Model.addConstr(SP::gamma2[n][t][i] == 0);
				//Model.addConstr(SP::beta[n][t][i] == SP::dual_val_beta[n][t][i]);
				/*if (EV::val_Xop[n][i] == 0)
				{
					Model.addConstr(SP::beta[n][t][i] == 0);
				}*/
				if (!(i == 5 || i == 8 || i == 9 || i == 10))
				{
					//Model.addConstr(SP::beta[n][t][i] == 0);
				}

			}
		}
	}

	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			double fge = 0;
			for (int n : Gnodes[k].adjE)
			{
				//Model.addConstr(SP::rho[k][n][tau] == 0); // commented
			}
		}
	}
	//Model.addConstr(SP::tau == 0);
	//Model.addConstr(SP::omega == 0);
#pragma endregion


#pragma region Objective Function
	GRBLinExpr ex0(0);

	//for (int i = 0; i < nPlt; i++)
	//{
	//	if (Plants[i].is_exis != 1) { continue; } // only applies to existing plants
	//	ex0 += Plants[i].max_yearly_gen * SP::alpha[i];
	//}

	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			double exp_prod = 0;
			for (int i = 0; i < nPlt; i++)
			{
				exp_prod += EV::val_prod[n][t][i];
			}
			double dem = Enodes[n].demand[Te[t]];
			double str = 0;
			for (int r = 0; r < neSt; r++)
			{
				str += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r];
			}
			ex0 += (dem - str - exp_prod - EV::val_curtE[n][t]) * SP::theta1[n][t];
			ex0 -= (dem - str - exp_prod - EV::val_curtE[n][t]) * SP::theta2[n][t];
			ex0 -= pi * SP::eta1[n][t];
			ex0 -= pi * SP::eta2[n][t];
		}
	}

	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (Branches[b].is_exist == 1)
			{
				ex0 -= Branches[b].maxFlow * SP::delta11[b][t];
				ex0 -= Branches[b].maxFlow * SP::delta12[b][t];
			}
			else
			{
				double Big_M = std::abs(Branches[b].suscep * 24);
				ex0 -= Branches[b].maxFlow * EV::val_Ze[b] * SP::delta21[b][t];
				ex0 -= Branches[b].maxFlow * EV::val_Ze[b] * SP::delta22[b][t];
				ex0 -= Big_M * (1 - EV::val_Ze[b]) * SP::zeta11[b][t];
				ex0 -= Big_M * (1 - EV::val_Ze[b]) * SP::zeta12[b][t];
			}
		}
	}



	//ex0 = 11;
	//ex0 = ex0 + SP::tau;

	Model.setObjective(ex0, GRB_MAXIMIZE);
	//Model.addConstr(ex0 == 9.68247e+08);
	//Model.addConstr(ex0 <= 0);
#pragma endregion


#pragma region Constraints for prod variable
	/*for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if ((Plants[i].type == "hydro-new" || Plants[i].type == "nuclear-new" ||
					Plants[i].type == "dfo" || Plants[i].type == "coal" ||
					Plants[i].type == "wind_offshore"))
				{
					continue;
				}


				GRBLinExpr ex_gam(0);
				GRBLinExpr ex_alpha(0);
				ex_alpha += 0;
				if (Plants[i].is_exis == 1)
				{
					ex_alpha += SP::alpha[i];
				}

				if (t == 0)
				{
					ex_gam += -SP::gamma1[n][t + 1][i] + SP::gamma2[n][t + 1][i];
				}
				else if (t == Te.size() - 1)
				{
					ex_gam += SP::gamma1[n][t][i] - SP::gamma2[n][t][i];
				}
				else
				{
					ex_gam += SP::gamma1[n][t][i] - SP::gamma1[n][t + 1][i] -
						SP::gamma2[n][t][i] + SP::gamma2[n][t + 1][i];
				}

				if (Plants[i].type == "nuclear")
				{
					Model.addConstr(time_weight[t] * ex_alpha + SP::beta[n][t][i] + SP::theta[n][t]
						<= time_weight[t] * (Plants[i].var_cost + Params::nuclear_price * Plants[i].heat_rate));
				}
				else if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
				{
					GRBLinExpr ex_adj(0);
					for (int k = 0; k < nGnode; k++)
					{
						for (int n2 : Gnodes[k].adjE)
						{
							if (n2 == n)
							{
								int tau = t / 24;
								ex_adj += Plants[i].heat_rate * SP::rho[k][n2][tau];
							}
						}
					}
					Model.addConstr(time_weight[t] * ex_alpha + SP::beta[n][t][i] +
						ex_gam + SP::theta[n][t] + ex_adj +
						time_weight[t] * Plants[i].emis_rate * Plants[i].heat_rate * SP::tau
						<= time_weight[t] * (Plants[i].var_cost));
				}
				else if (Plants[i].type == "solar" || Plants[i].type == "solar-UPV" ||
					Plants[i].type == "wind" || Plants[i].type == "wind-new" ||
					Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new")
				{
					Model.addConstr(time_weight[t] * ex_alpha + SP::beta[n][t][i]
						+ SP::theta[n][t] + SP::pi[n][t][i] + time_weight[t] * SP::omega
						<= time_weight[t] * (Plants[i].var_cost));
				}
				else if (Plants[i].type == "hydro" || Plants[i].type == "hydro-new")
				{
					Model.addConstr(time_weight[t] * ex_alpha + SP::beta[n][t][i]
						+ SP::theta[n][t] + time_weight[t] * SP::omega
						<= time_weight[t] * (Plants[i].var_cost));
				}
				else
				{
					Model.addConstr(time_weight[t] * ex_alpha + SP::beta[n][t][i]
						+ SP::theta[n][t]
						<= time_weight[t] * (Plants[i].var_cost));
				}
			}
		}
	}*/
#pragma endregion

#pragma region Constraints for flowE variable (flowE is a free variable)
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			GRBLinExpr ex_fe(0);
			if (Branches[b].is_exist == 1)
			{
				ex_fe -= SP::delta11[b][t];
				ex_fe += SP::delta12[b][t];
				ex_fe -= SP::zeta22[b][t];
				ex_fe += SP::zeta21[b][t];
			}
			else
			{
				ex_fe -= SP::delta21[b][t];
				ex_fe += SP::delta22[b][t];
				ex_fe -= SP::zeta11[b][t];
				ex_fe += SP::zeta12[b][t];
			}
			/*if (Branches[b].from_bus > Branches[b].to_bus)
			{
				ex_fe -= SP::theta[Branches[b].from_bus][t];
			}
			else
			{
				ex_fe += SP::theta[Branches[b].from_bus][t];
			}*/
			ex_fe += (SP::theta1[Branches[b].from_bus][t] - SP::theta1[Branches[b].to_bus][t]);
			ex_fe -= (SP::theta2[Branches[b].from_bus][t] - SP::theta2[Branches[b].to_bus][t]);
			//std::cout << "from bus: " << Branches[b].from_bus << " \t to bus: " << Branches[b].to_bus << endl;
			Model.addConstr(ex_fe == 0);
			//Model.addConstr(-ex_fe <= 1e-3);
		}
	}
#pragma endregion


#pragma region Constratints for phase angle, load shedding, 
	//// load shedding
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		Model.addConstr(SP::theta1[n][t] + SP::phi[n][t] <= time_weight[t] * E_curt_cost);
	//	}
	//}

	// phase angle (theta)
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			GRBLinExpr exp_line(0);

			for (int m : Enodes[n].adj_buses)
			{
				// instead of defining a dictionray, one could traverse the entire branches...
				int key1 = n * 200 + m;
				// check if this key exist, if not, the other order exisit
				if (Le.count(key1) == 0)
				{
					key1 = m * 200 + n;
				}
				//each Le can contain multiple lines
				for (int l2 : Le[key1])
				{
					if (Branches[l2].is_exist == 1)
					{
						if (n > m)
						{
							exp_line -= Branches[l2].suscep * (SP::zeta21[l2][t] - SP::zeta22[l2][t]);
						}
						else
						{
							exp_line += Branches[l2].suscep * (SP::zeta21[l2][t] - SP::zeta22[l2][t]);

						}
					}
					else
					{
						if (n > m)
						{
							exp_line += Branches[l2].suscep * (SP::zeta11[l2][t] - SP::zeta12[l2][t]);
						}
						else
						{
							exp_line -= Branches[l2].suscep * (SP::zeta11[l2][t] - SP::zeta12[l2][t]);
						}
					}
				}
			}

			if (t == 0)
			{
				Model.addConstr(exp_line - SP::eta1[n][t] + SP::eta2[n][t] + SP::eta3[t] == 0);
			}
			else
			{
				Model.addConstr(exp_line - SP::eta1[n][t] + SP::eta2[n][t] == 0);
			}
		}
	}
#pragma endregion
	Model.set(GRB_IntParam_OutputFlag, 0);
	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	//Model.set(GRB_IntParam_DualReductions, 0);
	Model.set(GRB_IntParam_InfUnbdInfo, 1);// if additional info for infeasible/unbounded models is desired

	Model.set(GRB_IntParam_Presolve, 0); // turn off presolve to be able to get extreme rays/farkas duals
	Model.set(GRB_IntParam_Method, 1); // 1: dual simplex, 0:simplex

	Model.optimize();
	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		//std::cout << "Failed to optimize IM!!!" << endl;
		//std::cout << Model.get(GRB_IntAttr_Status) << endl;;

		sp_new.dual_val_delta11 = new double* [nBr];
		sp_new.dual_val_delta12 = new double* [nBr];
		sp_new.dual_val_delta21 = new double* [nBr];
		sp_new.dual_val_delta22 = new double* [nBr];
		sp_new.dual_val_theta1 = new double* [nEnode];
		sp_new.dual_val_theta2 = new double* [nEnode];
		sp_new.dual_val_zeta11 = new double* [nBr];
		sp_new.dual_val_zeta12 = new double* [nBr];
		sp_new.dual_val_zeta21 = new double* [nBr];
		sp_new.dual_val_zeta22 = new double* [nBr];
		sp_new.dual_val_eta1 = new double* [nEnode];
		sp_new.dual_val_eta2 = new double* [nEnode];
		double fd = 0;
		for (int n = 0; n < nEnode; n++)
		{
			sp_new.dual_val_theta1[n] = new double[Te.size()]();
			sp_new.dual_val_theta2[n] = new double[Te.size()]();
			sp_new.dual_val_eta1[n] = new double[Te.size()]();
			sp_new.dual_val_eta2[n] = new double[Te.size()]();
			double ave_price = 0;
			for (int t = 0; t < Te.size(); t++)
			{
				sp_new.dual_val_theta1[n][t] = (std::max(0.0, SP::theta1[n][t].get(GRB_DoubleAttr_UnbdRay)));
				sp_new.dual_val_theta2[n][t] = (std::max(0.0, SP::theta2[n][t].get(GRB_DoubleAttr_UnbdRay)));
				//cout << "theta1: " << sp_new.dual_val_theta1[n][t] << " theta2: " << sp_new.dual_val_theta2[n][t] << endl;
				fd += std::abs(sp_new.dual_val_theta1[n][t]) + std::abs(sp_new.dual_val_theta2[n][t]);
				if (t > 0)
				{
					sp_new.dual_val_eta1[n][t] = std::max(0.0, SP::eta1[n][t].get(GRB_DoubleAttr_UnbdRay));
					fd += std::abs(sp_new.dual_val_eta1[n][t]);
					sp_new.dual_val_eta2[n][t] = std::max(0.0, SP::eta2[n][t].get(GRB_DoubleAttr_UnbdRay));
					fd += std::abs(sp_new.dual_val_eta2[n][t]);
					//cout << "eta1: " << sp_new.dual_val_eta1[n][t] << " eta2: " << sp_new.dual_val_eta2[n][t] << endl;

				}
			}
		}

		for (int b = 0; b < nBr; b++)
		{
			sp_new.dual_val_delta11[b] = new double[Te.size()]();
			sp_new.dual_val_delta12[b] = new double[Te.size()]();
			sp_new.dual_val_delta21[b] = new double[Te.size()]();
			sp_new.dual_val_delta22[b] = new double[Te.size()]();
			sp_new.dual_val_zeta11[b] = new double[Te.size()]();
			sp_new.dual_val_zeta12[b] = new double[Te.size()]();
			sp_new.dual_val_zeta21[b] = new double[Te.size()]();
			sp_new.dual_val_zeta22[b] = new double[Te.size()]();
			for (int t = 0; t < Te.size(); t++)
			{
				if (Branches[b].is_exist == 1)
				{
					sp_new.dual_val_delta11[b][t] = std::max(0.0, SP::delta11[b][t].get(GRB_DoubleAttr_UnbdRay));
					sp_new.dual_val_delta12[b][t] = std::max(0.0, SP::delta12[b][t].get(GRB_DoubleAttr_UnbdRay));
					sp_new.dual_val_zeta21[b][t] = std::max(0.0, SP::zeta21[b][t].get(GRB_DoubleAttr_UnbdRay));
					sp_new.dual_val_zeta22[b][t] = std::max(0.0, SP::zeta22[b][t].get(GRB_DoubleAttr_UnbdRay));
				}
				else
				{
					sp_new.dual_val_delta21[b][t] = std::max(0.0, SP::delta21[b][t].get(GRB_DoubleAttr_UnbdRay));
					sp_new.dual_val_delta22[b][t] = std::max(0.0, SP::delta22[b][t].get(GRB_DoubleAttr_UnbdRay));


					sp_new.dual_val_zeta11[b][t] = std::max(0.0, SP::zeta11[b][t].get(GRB_DoubleAttr_UnbdRay));
					sp_new.dual_val_zeta12[b][t] = std::max(0.0, SP::zeta12[b][t].get(GRB_DoubleAttr_UnbdRay));
				}
				//	cout << "existing: " << Branches[b].is_exist << " \t delta11: " << sp_new.dual_val_delta11[b][t] << " delta12: " << sp_new.dual_val_delta12[b][t] << " delta21: " << sp_new.dual_val_delta21[b][t] << " delta22: " << sp_new.dual_val_delta22[b][t] << endl;
					//cout << "zeta11: " << sp_new.dual_val_zeta11[b][t] << " zeta12: " << sp_new.dual_val_zeta12[b][t] << " zeta21: " << sp_new.dual_val_zeta21[b][t] << " zeta22: " << sp_new.dual_val_zeta22[b][t] << endl;

					//fd += SP::dual_val_delta11[b][t] + SP::dual_val_delta12[b][t] + SP::dual_val_delta21[b][t] + SP::dual_val_delta22[b][t];

			}
		}
		Cuts.push_back(sp_new);


		// print the values
		/*for (int n = 0; n < nEnode; n++)
		{
			for (int t = 0; t < Te.size(); t++)
			{
				if (sp_new.dual_val_theta1[n][t] > 0)
				{
					std::cout << "val_theta1[" << n << "][" << t << "] = " << sp_new.dual_val_theta1[n][t] << endl;
				}
				if (sp_new.dual_val_theta2[n][t] > 0)
				{
					std::cout << "val_theta2[" << n << "][" << t << "] = " << sp_new.dual_val_theta2[n][t] << endl;
				}

				if (sp_new.dual_val_eta1[n][t] > 0)
				{
					std::cout << "val_eta1[" << n << "][" << t << "] = " << sp_new.dual_val_eta1[n][t] << endl;
				}
				if (sp_new.dual_val_eta2[n][t] > 0)
				{
					std::cout << "val_eta2[" << n << "][" << t << "] = " << sp_new.dual_val_eta2[n][t] << endl;
				}
			}
		}

		for (int b = 0; b < nBr; b++)
		{
			for (int t = 0; t < Te.size(); t++)
			{
				if (sp_new.dual_val_delta11[b][t] > 0)
				{
					std::cout << "val_delta11[" << b << "][" << t << "] = " << sp_new.dual_val_delta11[b][t] << endl;
				}
				if (sp_new.dual_val_delta12[b][t] > 0)
				{
					std::cout << "val_delta12[" << b << "][" << t << "] = " << sp_new.dual_val_delta12[b][t] << endl;
				}
				if (sp_new.dual_val_delta21[b][t] > 0)
				{
					std::cout << "val_delta21[" << b << "][" << t << "] = " << sp_new.dual_val_delta21[b][t] << endl;
				}
				if (sp_new.dual_val_delta22[b][t] > 0)
				{
					std::cout << "val_delta22[" << b << "][" << t << "] = " << sp_new.dual_val_delta22[b][t] << endl;
				}

				if (sp_new.dual_val_zeta11[b][t] > 0)
				{
					std::cout << "val_zeta11[" << b << "][" << t << "] = " << sp_new.dual_val_zeta11[b][t] << endl;
				}
				if (sp_new.dual_val_zeta12[b][t] > 0)
				{
					std::cout << "val_zeta12[" << b << "][" << t << "] = " << sp_new.dual_val_zeta12[b][t] << endl;
				}
				if (sp_new.dual_val_zeta21[b][t] > 0)
				{
					std::cout << "val_zeta21[" << b << "][" << t << "] = " << sp_new.dual_val_zeta21[b][t] << endl;
				}
				if (sp_new.dual_val_zeta22[b][t] > 0)
				{
					std::cout << "val_zeta22[" << b << "][" << t << "] = " << sp_new.dual_val_zeta22[b][t] << endl;
				}
			}
		}*/

		return 0;

	}
	else
	{
		double obj_val = Model.get(GRB_DoubleAttr_ObjVal);

		//double gap = Model.get(GRB_DoubleAttr_MIPGap);

		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		std::cout << "\n\n\n\n Dual Elapsed time: " << Elapsed << endl;
		std::cout << "\t Dual Obj Value:" << obj_val << endl;
		std::cout << " Dual Status:" << Model.get(GRB_IntAttr_Status) << endl;
		return -1;
	}
	/*if (~Model.get(GRB_IntAttr_Status)==GRB_OPTIMAL)
	{
		std::cout << "Failed to optimize IM!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
	}*/




#pragma region print the variables

	//// alpha
	//for (int i = 0; i < nPlt; i++)
	//{
	//	double alpha = SP::alpha[i].get(GRB_DoubleAttr_X);
	//	if (alpha < -0.001 || alpha>0.001)
	//	{
	//		std::cout << "alpha[" << i << "] = " << alpha << endl;
	//	}
	//}

	//// beta
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (t > 5)
	//		{
	//			break;
	//		}
	//		for (int i = 0; i < nPlt; i++)
	//		{
	//			double beta = SP::beta[n][t][i].get(GRB_DoubleAttr_X);
	//			if (beta < -0.001 || beta>0.001)
	//			{
	//				std::cout << "beta[" << n << "][" << t << "][" << i << "] = " << beta << endl;
	//			}
	//		}
	//	}
	//}

	////gamma
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (t > 10)
	//		{
	//			break;
	//		}
	//		for (int i = 0; i < nPlt; i++)
	//		{
	//			double gamma = SP::gamma1[n][t][i].get(GRB_DoubleAttr_X);
	//			if (gamma < -0.001 || gamma>0.001)
	//			{
	//				std::cout << "gamma1[" << n << "][" << t << "][" << i << "] = " << gamma << endl;
	//			}

	//			gamma = SP::gamma2[n][t][i].get(GRB_DoubleAttr_X);
	//			if (gamma < -0.001 || gamma>0.001)
	//			{
	//				std::cout << "gamma2[" << n << "][" << t << "][" << i << "] = " << gamma << endl;
	//			}
	//		}
	//	}
	//}


	//theta
	/*for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (t > 10)
			{
				break;
			}
			double theta = SP::theta1[n][t].get(GRB_DoubleAttr_X);
			if (theta > 0.001 || theta < -0.001)
			{
				std::cout << "theta1[" << n << "][" << t << "] = " << theta << endl;
			}
			theta = SP::theta2[n][t].get(GRB_DoubleAttr_X);
			if (theta > 0.001 || theta < -0.001)
			{
				std::cout << "theta2[" << n << "][" << t << "] = " << theta << endl;
			}
		}
	}*/

	////// alpha from primal subproblem
	////for (int i = 0; i < nPlt; i++)
	////{
	////	double alpha = EV::dual_val_alpha[i];
	////	if (alpha < -0.001 || alpha>0.001)
	////	{
	////		std::cout << "PSP_alpha[" << i << "] = " << alpha << endl;
	////	}
	////}

	////// beta from primal subproblem
	////for (int n = 0; n < nEnode; n++)
	////{
	////	for (int t = 0; t < Te.size(); t++)
	////	{
	////		if (t > 10)
	////		{
	////			break;
	////		}
	////		for (int i = 0; i < nPlt; i++)
	////		{
	////			double beta = EV::dual_val_beta[n][t][i];
	////			if (!(i == 5 || i == 8 || i == 9 || i == 10)) { continue; }
	////			if (beta < -0.001 || beta>0.001)
	////			{
	////				std::cout << "PSP_beta[" << n << "][" << t << "][" << i << "] = " << beta << endl;
	////			}
	////		}
	////	}
	////}


	//// theta from primal subproblem
	///*for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (t > 10)
	//		{
	//			break;
	//		}
	//		double theta = SP::dual_val_theta[n][t];
	//		if (theta > 0.001 || theta < -0.001)
	//		{
	//			std::cout << "PSP_theta[" << n << "][" << t << "] = " << theta << endl;
	//		}
	//	}
	//}*/

	////// phi from primal subproblem
	////for (int n = 0; n < nEnode; n++)
	////{
	////	for (int t = 0; t < Te.size(); t++)
	////	{
	////		if (t > 10)
	////		{
	////			break;
	////		}
	////		double phi = EV::dual_val_phi[n][t];
	////		//if (theta > 0.001 || theta < -0.001)
	////		{
	////			std::cout << "PSP_phi[" << n << "][" << t << "] = " << phi << endl;
	////		}
	////	}
	////}



	////phi
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (t > 10)
	//		{
	//			//break;
	//		}
	//		double phi = SP::phi[n][t].get(GRB_DoubleAttr_X);
	//		if (phi > 0.001 || phi < -0.001)
	//		{
	//			std::cout << "phi[" << n << "][" << t << "] = " << phi << endl;
	//		}
	//	}
	//}

	////pi
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (t > 10)
	//		{
	//			break;
	//		}
	//		for (int i = 0; i < nPlt; i++)
	//		{
	//			double pi = SP::pi[n][t][i].get(GRB_DoubleAttr_X);
	//			if (pi < -0.001 || pi>0.001)
	//			{
	//				std::cout << "pi[" << n << "][" << t << "][" << i << "] = " << pi << endl;
	//			}
	//		}
	//	}
	//}

	//// delta
	//for (int b = 0; b < nBr; b++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (t > 10)
	//		{
	//			//break;
	//		}
	//		double dz = SP::delta11[b][t].get(GRB_DoubleAttr_X);
	//		if (dz < -0.001 || dz>0.001)
	//		{
	//			std::cout << "delta11[" << b << "][" << t << "] = " << dz << endl;
	//		}
	//		dz = SP::delta12[b][t].get(GRB_DoubleAttr_X);
	//		if (dz < -0.001 || dz>0.001)
	//		{
	//			std::cout << "delta12[" << b << "][" << t << "] = " << dz << endl;
	//		}

	//		dz = SP::delta21[b][t].get(GRB_DoubleAttr_X);
	//		if (dz < -0.001 || dz>0.001)
	//		{
	//			std::cout << "delta21[" << b << "][" << t << "] = " << dz << endl;
	//		}

	//		dz = SP::delta22[b][t].get(GRB_DoubleAttr_X);
	//		if (dz < -0.001 || dz>0.001)
	//		{
	//			std::cout << "delta22[" << b << "][" << t << "] = " << dz << endl;
	//		}
	//	}
	//}

	//////zeta
	//for (int b = 0; b < nBr; b++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (t > 10)
	//		{
	//			//break;
	//		}
	//		double dz = SP::zeta11[b][t].get(GRB_DoubleAttr_X);
	//		if (dz < -0.001 || dz>0.001)
	//		{
	//			std::cout << "zeta11[" << b << "][" << t << "] = " << dz << endl;
	//		}
	//		dz = SP::zeta12[b][t].get(GRB_DoubleAttr_X);
	//		if (dz < -0.001 || dz>0.001)
	//		{
	//			std::cout << "zeta12[" << b << "][" << t << "] = " << dz << endl;
	//		}
	//		dz = SP::zeta21[b][t].get(GRB_DoubleAttr_X);
	//		if (dz < -0.001 || dz>0.001)
	//		{
	//			std::cout << "zeta21[" << b << "][" << t << "] = " << dz << endl;
	//		}
	//		dz = SP::zeta22[b][t].get(GRB_DoubleAttr_X);
	//		if (dz < -0.001 || dz>0.001)
	//		{
	//			std::cout << "zeta22[" << b << "][" << t << "] = " << dz << endl;
	//		}
	//	}
	//}


	////eta
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (t > 10)
	//		{
	//			//break;
	//		}
	//		double eta = SP::eta1[n][t].get(GRB_DoubleAttr_X);
	//		if (eta < -0.001 || eta>0.001)
	//		{
	//			std::cout << "eta1[" << n << "][" << t << "] = " << eta << endl;
	//		}

	//		eta = SP::eta2[n][t].get(GRB_DoubleAttr_X);
	//		if (eta < -0.001 || eta>0.001)
	//		{
	//			std::cout << "eta1[" << n << "][" << t << "] = " << eta << endl;
	//		}

	//		eta = SP::eta3[t].get(GRB_DoubleAttr_X);
	//		if (eta < -0.001 || eta>0.001)
	//		{
	//			std::cout << "eta1[" << n << "][" << t << "] = " << eta << endl;
	//		}
	//	}
	//}

	////rho

	//for (int k = 0; k < nGnode; k++)
	//{
	//	for (int n : Gnodes[k].adjE)
	//	{
	//		for (int tau = 0; tau < Tg.size(); tau++)
	//		{
	//			double rho = SP::rho[k][n][tau].get(GRB_DoubleAttr_X);
	//			if (rho < -0.001 || rho>0.001)
	//			{
	//				std::cout << "rho[" << k << "][" << n << "][" << tau << "] = " << rho << endl;
	//			}
	//		}
	//	}
	//}

	// print omega
	//std::cout << "omega = " << SP::omega.get(GRB_DoubleAttr_X) << endl;
	//std::cout << "tau = " << SP::tau.get(GRB_DoubleAttr_X) << endl;
#pragma endregion


#pragma region verify the obj value


	//double obj_theta = 0; double obj_other = 0;
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{

	//		double dem = Enodes[n].demand[Te[t]];
	//		double str = 0;
	//		for (int r = 0; r < neSt; r++)
	//		{
	//			str += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r];
	//		}
	//		//std::cout << "dem-str: " << dem - str << endl;
	//		obj_theta += (dem - str) * SP::theta1[n][t].get(GRB_DoubleAttr_X); // commented
	//		obj_other += pi * SP::eta1[n][t].get(GRB_DoubleAttr_X) + pi * SP::eta2[n][t].get(GRB_DoubleAttr_X);
	//	}
	//}

	//for (int b = 0; b < nBr; b++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (b == 27 && t == 16)
	//		{
	//			std::cout << EV::val_Ze[b] << endl;
	//		}
	//		if (Branches[b].is_exist == 1)
	//		{
	//			obj_other += Branches[b].maxFlow * (SP::delta11[b][t].get(GRB_DoubleAttr_X) + SP::delta12[b][t].get(GRB_DoubleAttr_X));
	//		}
	//		else
	//		{
	//			obj_other += Branches[b].maxFlow * EV::val_Ze[b] * (SP::delta21[b][t].get(GRB_DoubleAttr_X) + SP::delta22[b][t].get(GRB_DoubleAttr_X));
	//			obj_other += std::abs(Branches[b].suscep) * 24 * (1 - EV::val_Ze[b]) * (SP::zeta11[b][t].get(GRB_DoubleAttr_X) + SP::zeta12[b][t].get(GRB_DoubleAttr_X));
	//		}
	//	}
	//}




	//std::cout << "dual obj_theta value: " << obj_theta << endl;
	//std::cout << "dual obj_other value: " << obj_other << endl;
	//std::cout << "total dual obj value: " << obj_theta + obj_other << endl;
#pragma endregion

}


double MP_init_heuristic(GRBEnv* env)
{
	// replace UC commitment constraints and variables
	// relax integrality of Xop, Xest, Xdec
	// solve and get a lower bound soluton

#pragma region Fetch Data
	std::map<string, int> sym2pltType = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7},
		{"CT",8},{"CC",9},{"CC-CCS",10},{"solar-UPV",11},{"wind-new",12},
		{"wind-offshore-new",13},{"hydro-new",14},{"nuclear-new",15} };
	std::map<int, string> pltType2sym = { {0,"ng"},{1,"dfo"},
{2,"solar"},{3,"wind"},{4,"wind_offshore"},{5,"hydro"},{6,"coal"},{7,"nuclear"} };
	vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Params::Enodes.size();
	int nPlt = (int)Params::Plants.size();
	int nBr = (int)Params::Branches.size();
	int neSt = (int)Params::Estorage.size();
	vector<int> Tg = Params::Tg;
	vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	double pi = 3.141592;
	int nGnode = (int)Params::Gnodes.size();
	double WACC = Params::WACC;
	int trans_unit_cost = Params::trans_unit_cost;
	int trans_line_lifespan = Params::trans_line_lifespan;
	double NG_price = Params::NG_price;
	double dfo_pric = Params::dfo_pric;
	double coal_price = Params::coal_price;
	double nuclear_price = Params::nuclear_price;
	double E_curt_cost = Params::E_curt_cost;
	double G_curt_cost = Params::G_curt_cost;
	double pipe_per_mile = Params::pipe_per_mile;
	int pipe_lifespan = Params::pipe_lifespan;
	int battery_lifetime = Params::battery_lifetime;
	map<int, vector<int>> Le = Params::Le;
	vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion
	Setting::UC_active = false; // only applies to the full problem
	Setting::relax_int_vars = true; // int vars (but not binary vars) in electricity network

	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_NGobj(0);
	GRBLinExpr exp_Eobj(0);
	Populate_EV_SP(Model);
	Elec_Module(Model, exp_Eobj);

	Populate_GV(Model);
	NG_Module(Model, exp_NGobj);


	Setting::UC_active = true; // only applies to the full problem
	Setting::relax_int_vars = false; // int vars (but not binary vars) in electricity network

	// coupling constraints
	GRBLinExpr ex_xi(0);
	GRBLinExpr ex_NG_emis(0);
	GRBLinExpr ex_E_emis(0);
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);


	Model.addConstr(CV::xi == ex_xi);

	if (Setting::Case == 3)
	{// the original model
		Model.addConstr(ex_E_emis + ex_NG_emis <= Setting::Emis_redu_goal * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}

	Model.setObjective(exp_Eobj + exp_NGobj, GRB_MINIMIZE);

#pragma region Solve the model
	//Model.set(GRB_IntParam_OutputFlag, 0);
	Model.set(GRB_DoubleParam_TimeLimit, 600);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	//Model.set(GRB_IntParam_DualReductions, 0);

	Model.optimize();
	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to LB problem !!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
		return -1;
	}
	double obj_val = Model.get(GRB_DoubleAttr_ObjBound);
	//double obj_val = Model.get(GRB_DoubleAttr_ObjVal);
	std::cout << "LB obj: " << obj_val;

#pragma endregion

	return obj_val;
}

void MP_init_heuristic2(GRBEnv* env)
{
	// this is what I tried before but didn't seem to work good. 
	// this method generates a feasiblt solution to be fed to the MP in the first iteration
	// as a warm start solution
	//auto start = chrono::high_resolution_clock::now();
#pragma region Fetch Data
	std::map<string, int> sym2pltType = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7},
		{"CT",8},{"CC",9},{"CC-CCS",10},{"solar-UPV",11},{"wind-new",12},
		{"wind-offshore-new",13},{"hydro-new",14},{"nuclear-new",15} };
	std::map<int, string> pltType2sym = { {0,"ng"},{1,"dfo"},
{2,"solar"},{3,"wind"},{4,"wind_offshore"},{5,"hydro"},{6,"coal"},{7,"nuclear"} };
	vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Params::Enodes.size();
	int nPlt = (int)Params::Plants.size();
	int nBr = (int)Params::Branches.size();
	int neSt = (int)Params::Estorage.size();
	vector<int> Tg = Params::Tg;
	vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	double pi = 3.141592;
	int nGnode = (int)Params::Gnodes.size();
	double WACC = Params::WACC;
	int trans_unit_cost = Params::trans_unit_cost;
	int trans_line_lifespan = Params::trans_line_lifespan;
	double NG_price = Params::NG_price;
	double dfo_pric = Params::dfo_pric;
	double coal_price = Params::coal_price;
	double nuclear_price = Params::nuclear_price;
	double E_curt_cost = Params::E_curt_cost;
	double G_curt_cost = Params::G_curt_cost;
	double pipe_per_mile = Params::pipe_per_mile;
	int pipe_lifespan = Params::pipe_lifespan;
	int battery_lifetime = Params::battery_lifetime;
	map<int, vector<int>> Le = Params::Le;
	vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion


#pragma region Firts phase: ignore UC vars and network cons.->solve->get Xop vars

#pragma region Initialization: NG module, define vars, etc.
	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_NGobj(0);
	GRBLinExpr exp_Eobj(0);

	Populate_EV_SP(Model);
	Populate_GV(Model);
	NG_Module(Model, exp_NGobj);

	// coupling constraints
	GRBLinExpr ex_xi(0);
	GRBLinExpr ex_NG_emis(0);
	GRBLinExpr ex_E_emis(0);
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);

	//if (Setting::Case == 1)
	//{
	//	Model.addConstr(CV::xi == 100); // just to populate CV::xi
	//} else if to next 

	if (Setting::is_xi_given)
	{
		//std::cout << "\n\n if xi given" << endl;
		Model.addConstr(ex_xi == Setting::xi_val * Setting::PGC);
		Model.addConstr(CV::xi == ex_xi);
	}
	else
	{
		//std::cout << "else xi given" << endl;
		Model.addConstr(CV::xi == ex_xi);
	}

	/*if (Setting::Case == 1)
	{
		Model.addConstr(100 == CV::E_emis);
		Model.addConstr(100 == CV::NG_emis);
	}*/
	if (Setting::Case == 2)
	{
		Model.addConstr(ex_E_emis <=(1- Setting::Emis_redu_goal) * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}
	if (Setting::Case == 3)
	{// the original model
		Model.addConstr(ex_E_emis + ex_NG_emis <=(1- Setting::Emis_redu_goal) * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}
#pragma endregion

#pragma region Set some variables
	// 1) existing types can not be established because there are new equivalent types
	// 2) new types cannot be decommissioned
	// 3) no new wind-offshore on some location (not yet implemented)
	// 4) no new nuclear
	// 5) turn off df0/coal/wind-offshore. The data is not correct on wind-offshore, there is no existing wind-offshore in the region
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			//// do not allow establishment off-shore wind and hydro plants
			//Model.addConstr(EV::Xdec[n][5] == 0);
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore")
			{
				Model.addConstr(EV::Xop[n][i] == 0);
			}

			if (Plants[i].type == "hydro-new")
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}

			if (Plants[i].type == "wind-offshore-new")
			{
				if (n == 0 || n == 2)// except for Massachusetts and Connecticus
				{
					continue;
				}
				Model.addConstr(EV::Xest[n][i] == 0);
			}

			/*if (Plants[i].type == "nuclear-new")
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}*/

			if (Plants[i].is_exis == 1)
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}
			else
			{
				Model.addConstr(EV::Xdec[n][i] == 0);
			}


		}
	}
#pragma endregion


#pragma region Objective Function
	GRBLinExpr ex_est(0);
	GRBLinExpr ex_decom(0);
	GRBLinExpr ex_fix(0);
	GRBLinExpr ex_startup(0);
	GRBLinExpr ex_var(0);
	GRBLinExpr ex_thermal_fuel(0);
	GRBLinExpr ex_emis(0);
	GRBLinExpr ex_shedd(0);
	GRBLinExpr ex_trans(0);
	GRBLinExpr ex_elec_str(0);


	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			// Investment/decommission cost of plants
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore") {
				continue;
			}
			double s1 = (double)std::pow(1.0 / (1 + WACC), Plants[i].lifetime);
			double capco = WACC / (1 - s1) * Plants[i].Reg_coeffs_per_state[n];
			ex_est += capco * Plants[i].capex * EV::Xest[n][i];
			ex_decom += Plants[i].decom_cost * EV::Xdec[n][i];
		}

		// fixed cost (annual, so no iteration over time)
		for (int i = 0; i < nPlt; i++)
		{
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore") {
				continue;
			}
			ex_fix += Plants[i].fix_cost * EV::Xop[n][i];
		}
		// var+fuel costs of plants
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
					Plants[i].type == "wind_offshore") {
					continue;
				}
				// var cost
				ex_var += time_weight[t] * Plants[i].var_cost * EV::prod[n][t][i];

				// fuel price to be updated later (dollar per thousand cubic feet=MMBTu)
				// NG fuel is handle in the NG network for case 2 and 3
				//if (Plants[i].type == "dfo")
				//{
				//	ex_thermal_fuel += time_weight[t] * dfo_pric * Plants[i].heat_rate * EV::prod[n][t][i];

				//}
				//else if (Plants[i].type == "coal")
				//{
				//	ex_thermal_fuel += time_weight[t] * coal_price * Plants[i].heat_rate * EV::prod[n][t][i];
				//}
				if (Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
				{
					ex_thermal_fuel += time_weight[t] * nuclear_price * Plants[i].heat_rate * EV::prod[n][t][i];
				}
				else if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
				{
					if (Setting::Case == 1)
					{
						ex_thermal_fuel += time_weight[t] * NG_price * Plants[i].heat_rate * EV::prod[n][t][i];
					}
				}

				//if (Plants[i].type == "ng" || Plants[i].type == "CT" ||
				//	Plants[i].type == "CC" || Plants[i].type == "CC-CCS" ||
				//	Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
				//{
				//	// startup cost
				//	ex_startup += time_weight[t] * Plants[i].startup_cost * EV::Xup[n][t][i];
				//}

			}

			// load curtailment cost
			ex_shedd += time_weight[t] * E_curt_cost * EV::curtE[n][t];
		}


		// storage cost
		for (int r = 0; r < neSt; r++)
		{
			double s1 = (double)std::pow(1.0 / (1 + WACC), battery_lifetime);
			double capco = WACC / (1 - s1);
			ex_elec_str += capco * Estorage[r].power_cost * EV::YeCD[n][r] + capco * Estorage[r].energy_cost * EV::YeLev[n][r];
			ex_elec_str += Estorage[r].pFOM * EV::YeCD[n][r] + Estorage[r].eFOM * EV::YeLev[n][r]; // fixed cost per kw per year
		}
	}

	// investment cost of transmission lines
	for (int b = 0; b < Branches.size(); b++)
	{
		double s1 = (double)std::pow(1.0 / (1 + WACC), trans_line_lifespan);
		double capco = WACC / (1 - s1);
		int fb = Branches[b].from_bus;
		int tb = Branches[b].to_bus;
		// find the index of tb in the adj_buses of the fb
		int tbi = std::find(Enodes[fb].adj_buses.begin(), Enodes[fb].adj_buses.end(), tb) - Enodes[fb].adj_buses.begin();
		ex_trans += capco * trans_unit_cost * Branches[b].maxFlow * Branches[b].length * EV::Ze[b];
	}
	exp_Eobj = ex_est + ex_decom + ex_fix + ex_startup + ex_emis + ex_var + ex_thermal_fuel + ex_shedd + ex_trans + ex_elec_str;

	Model.addConstr(EV::est_cost == ex_est);
	Model.addConstr(EV::est_trans_cost == ex_trans);
	Model.addConstr(EV::decom_cost == ex_decom);
	Model.addConstr(EV::fixed_cost == ex_fix);
	Model.addConstr(EV::var_cost == ex_var);
	Model.addConstr(EV::thermal_fuel_cost == ex_thermal_fuel);
	Model.addConstr(EV::dfo_coal_emis_cost == ex_emis);
	Model.addConstr(EV::shedding_cost == ex_shedd);
	Model.addConstr(EV::elec_storage_cost == ex_elec_str);
	Model.addConstr(EV::startup_cost == ex_startup);

	Model.addConstr(exp_Eobj == EV::e_system_cost);
	Model.setObjective(exp_Eobj + exp_NGobj, GRB_MINIMIZE);
#pragma endregion


#pragma region MP Electricity local constraints
	// C1: number of generation units at each node	
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			// is plant type i part of initial plants types of node n:
			bool isInd = std::find(Enodes[n].Init_plt_ind.begin(), Enodes[n].Init_plt_ind.end(), i) != Enodes[n].Init_plt_ind.end();
			if (isInd)
			{
				// find the index of plant in the set of Init_plants of the node
				string tp = pltType2sym[i];
				int ind1 = std::find(Enodes[n].Init_plt_types.begin(), Enodes[n].Init_plt_types.end(), tp) - Enodes[n].Init_plt_types.begin();
				Model.addConstr(EV::Xop[n][i] == Enodes[n].Init_plt_count[ind1] - EV::Xdec[n][i] + EV::Xest[n][i]);
			}
			else
			{
				Model.addConstr(EV::Xop[n][i] == -EV::Xdec[n][i] + EV::Xest[n][i]);
			}

			for (int t = 0; t < Te.size(); t++)
			{
				// (w/o UC version) keep this to avoid power generation from hydro-new, dfor, etc.
				Model.addConstr(EV::prod[n][t][i] <= Plants[i].Pmax * EV::Xop[n][i]);
				// (w/o UC version)
				if (t > 0)
				{
					//Model.addConstr(Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i] >= -EV::prod[n][t][i] + EV::prod[n][t - 1][i]);
					Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
					Model.addConstr(-EV::prod[n][t][i] + EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
				}
			}
		}
	}

	////C2, C3, C4, C5: UC,  production limit, ramping for thermal units (ng, CT, CC, CC-CCS, nuclear)	
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int i = 0; i < nPlt; i++)
	//	{
	//		for (int t = 0; t < Te.size(); t++)
	//		{
	//			if (Plants[i].type == "ng" || Plants[i].type == "CT" ||
	//				Plants[i].type == "CC" || Plants[i].type == "CC-CCS" ||
	//				Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
	//			{
	//				// UC 1
	//				if (t > 0)
	//				{
	//					Model.addConstr(EV::X[n][t][i] - EV::X[n][t - 1][i] == EV::Xup[n][t][i] - EV::Xdown[n][t][i]);
	//				}
	//				// UC 2
	//				Model.addConstr(EV::X[n][t][i] <= EV::Xop[n][i]);
	//				// output limit
	//				Model.addConstr(EV::prod[n][t][i] >= Plants[i].Pmin * Plants[i].Pmax * EV::X[n][t][i]);
	//				Model.addConstr(EV::prod[n][t][i] <= Plants[i].Pmax * EV::X[n][t][i]);
	//				//ramping
	//				if (t > 0)
	//				{
	//					double mx1 = std::max(Plants[i].Pmin, Plants[i].rampU);
	//					Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * (EV::X[n][t][i] - EV::Xup[n][t][i]) + mx1 * Plants[i].Pmax * EV::Xup[n][t][i]);
	//					Model.addConstr(-EV::prod[n][t][i] + EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * (EV::X[n][t][i] - EV::Xup[n][t][i]) + mx1 * Plants[i].Pmax * EV::Xup[n][t][i]);
	//				}
	//			}
	//		}
	//	}
	//}

	// C5, C6: flow limit for electricity
	//C7: power balance
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			GRBLinExpr exp_prod(0);
			for (int i = 0; i < nPlt; i++)
			{
				exp_prod += EV::prod[n][t][i];
			}
			GRBLinExpr exp_trans(0);
			for (int m : Enodes[n].adj_buses)
			{
				// instead of defining a dictionray, one could traverse the entire branches...
				int key1 = n * 200 + m;
				// check if this key exist, if not, the other order exisit
				if (Le.count(key1) == 0)
				{
					key1 = m * 200 + n;
				}
				//each Le can contain multiple lines
				for (int l2 : Le[key1])
				{
					if (n > m)
					{
						exp_trans -= EV::flowE[l2][t];
					}
					else
					{
						exp_trans += EV::flowE[l2][t];
					}

				}
			}
			GRBLinExpr ex_store(0);
			for (int r = 0; r < neSt; r++)
			{
				ex_store += EV::eSdis[n][t][r] - EV::eSch[n][t][r];
			}
			double dem = Enodes[n].demand[Te[t]];


			// ignore trans
			Model.addConstr(exp_prod + ex_store + EV::curtE[n][t] == dem);
			//Model.addConstr(exp_prod + exp_trans + ex_store + EV::curtE[n][t] == dem); //uncomment this line
		}
	}


	// C8: flow equation
	// C9: phase angle (theta) limits. already applied in the definition of the variable	
	// C10: VRE production profile
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "solar" || Plants[i].type == "solar-UPV")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::Xop[n][i]);
				}
				if (Plants[i].type == "wind" || Plants[i].type == "wind-new")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::Xop[n][i]);
				}
				if (Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][0] * Plants[i].Pmax * EV::Xop[n][i]);
				}
			}
		}
	}


	// C11: demand curtainlment constraint
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			double dem = Enodes[n].demand[Te[t]];
			//dem = dem/1.896;
			Model.addConstr(EV::curtE[n][t] <= dem);
		}
	}


	//C12: RPS constraints
	GRBLinExpr exp3(0);
	double total_demand = 0;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			total_demand += time_weight[t] * Enodes[n].demand[Te[t]];

			for (int i = 0; i < nPlt; i++)
			{
				//exp2 += time_weight[t] * Plants[i].emis_rate * EV::prod[n][t][i];
				if (Plants[i].type == "solar" || Plants[i].type == "wind" ||
					Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new" ||
					Plants[i].type == "hydro" || Plants[i].type == "solar-UPV" ||
					Plants[i].type == "wind-new" || Plants[i].type == "hydro-new")
				{
					exp3 += time_weight[t] * EV::prod[n][t][i]; // production from VRE
				}
			}
		}
	}
	Model.addConstr(exp3 >= Setting::RPS * total_demand);

	// C14,C15,C16 storage constraints
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int r = 0; r < neSt; r++)
			{
				/*if (t == Te.size() - 1)
				{
					Model.addConstr(EV::eSlev[n][t][r] <= EV::eSlev[n][0][r]);
				}*/
				if (t == 0)
				{
					Model.addConstr(EV::eSlev[n][t][r] ==
						Estorage[r].eff_ch * EV::eSch[n][t][r] -
						(EV::eSdis[n][t][r] / Estorage[r].eff_disCh));
				}
				else
				{
					Model.addConstr(EV::eSlev[n][t][r] == EV::eSlev[n][t - 1][r] +
						Estorage[r].eff_ch * EV::eSch[n][t][r] -
						(EV::eSdis[n][t][r] / Estorage[r].eff_disCh));
				}
				Model.addConstr(EV::YeCD[n][r] >= EV::eSdis[n][t][r]);
				Model.addConstr(EV::YeCD[n][r] >= EV::eSch[n][t][r]);
				Model.addConstr(EV::YeLev[n][r] >= EV::eSlev[n][t][r]);
			}
		}
	}

	// start and ending of storage should be the same in case of using rep. days
	for (int n = 0; n < nEnode; n++)
	{
		for (int r = 0; r < neSt; r++)
		{
			for (int k = 0; k < Tg.size(); k++)
			{
				int str = k * 24;
				int end = (k + 1) * 24 - 1;
				Model.addConstr(EV::eSlev[n][str][r] == EV::eSlev[n][end][r]);
			}
		}
	}
	// C17: if an storage established or not
	for (int n = 0; n < nEnode; n++)
	{
		for (int r = 0; r < neSt; r++)
		{
			Model.addConstr(EV::YeLev[n][r] <= 10e8 * EV::YeStr[n][r]);
			Model.addConstr(EV::YeCD[n][r] <= EV::YeLev[n][r]);
		}
	}
#pragma endregion

#pragma region Solve and get values
	Model.set(GRB_IntParam_OutputFlag, 0);
	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model.optimize();
	std::cout << "First phase obj: " << Model.get(GRB_DoubleAttr_ObjVal);
	EV::val_Xop = new double* [nEnode];
	EV::val_Xest = new double* [nEnode];
	EV::val_Xdec = new double* [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_Xop[n] = new double[nPlt]();
		EV::val_Xest[n] = new double[nPlt]();
		EV::val_Xdec[n] = new double[nPlt]();
		for (int i = 0; i < nPlt; i++)
		{
			EV::val_Xop[n][i] = std::round(EV::Xop[n][i].get(GRB_DoubleAttr_X));
			EV::val_Xest[n][i] = std::round(EV::Xest[n][i].get(GRB_DoubleAttr_X));
			EV::val_Xdec[n][i] = std::round(EV::Xdec[n][i].get(GRB_DoubleAttr_X));
		}
	}
#pragma endregion



#pragma endregion

#pragma region Second phase: fix Xop -> solve full problem -> get values of UC vars (X) and Xopt to feed MP

#pragma region Initialization: NG module, define vars, etc.
	//Model.terminate();
	Model.reset();
	//Model = new GRBModel();
	//GRBModel Model = GRBModel(env);
	exp_NGobj = 0;
	exp_Eobj = 0;

	Populate_EV_SP(Model);
	Populate_GV(Model);
	NG_Module(Model, exp_NGobj);

	// coupling constraints
	ex_xi = 0;
	ex_NG_emis = 0;
	ex_E_emis = 0;
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);

	//if (Setting::Case == 1)
	//{
	//	Model.addConstr(CV::xi == 100); // just to populate CV::xi
	//} else if to next 

	if (Setting::is_xi_given)
	{
		//std::cout << "\n\n if xi given" << endl;
		Model.addConstr(ex_xi == Setting::xi_val * Setting::PGC);
		Model.addConstr(CV::xi == ex_xi);
	}
	else
	{
		//std::cout << "else xi given" << endl;
		Model.addConstr(CV::xi == ex_xi);
	}

	/*if (Setting::Case == 1)
	{
		Model.addConstr(100 == CV::E_emis);
		Model.addConstr(100 == CV::NG_emis);
	}*/
	if (Setting::Case == 2)
	{
		Model.addConstr(ex_E_emis <=(1- Setting::Emis_redu_goal) * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}
	if (Setting::Case == 3)
	{// the original model
		Model.addConstr(ex_E_emis + ex_NG_emis <=(1- Setting::Emis_redu_goal) * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}
#pragma endregion

#pragma region Set some variables
	// 1) existing types can not be established because there are new equivalent types
	// 2) new types cannot be decommissioned
	// 3) no new wind-offshore on some location (not yet implemented)
	// 4) no new nuclear
	// 5) turn off df0/coal/wind-offshore. The data is not correct on wind-offshore, there is no existing wind-offshore in the region
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			//// do not allow establishment off-shore wind and hydro plants
			//Model.addConstr(EV::Xdec[n][5] == 0);
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore")
			{
				Model.addConstr(EV::Xop[n][i] == 0);
			}

			if (Plants[i].type == "hydro-new")
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}

			if (Plants[i].type == "wind-offshore-new")
			{
				if (n == 0 || n == 2)// except for Massachusetts and Connecticus
				{
					continue;
				}
				Model.addConstr(EV::Xest[n][i] == 0);
			}

			/*if (Plants[i].type == "nuclear-new")
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}*/

			if (Plants[i].is_exis == 1)
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}
			else
			{
				Model.addConstr(EV::Xdec[n][i] == 0);
			}


		}
	}

	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			Model.addConstr(EV::Xest[n][i] == EV::val_Xest[n][i]);
			Model.addConstr(EV::Xdec[n][i] == EV::val_Xdec[n][i]);
			Model.addConstr(EV::Xop[n][i] == EV::val_Xop[n][i]);
		}
	}

#pragma endregion


#pragma region Objective Function
	ex_est = 0;
	ex_decom = 0;
	ex_fix = 0;
	ex_startup = 0;
	ex_var = 0;
	ex_thermal_fuel = 0;
	ex_emis = 0;
	ex_shedd = 0;
	ex_trans = 0;
	ex_elec_str = 0;


	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			// Investment/decommission cost of plants
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore") {
				continue;
			}
			double s1 = (double)std::pow(1.0 / (1 + WACC), Plants[i].lifetime);
			double capco = WACC / (1 - s1) * Plants[i].Reg_coeffs_per_state[n];
			ex_est += capco * Plants[i].capex * EV::Xest[n][i];
			ex_decom += Plants[i].decom_cost * EV::Xdec[n][i];
		}

		// fixed cost (annual, so no iteration over time)
		for (int i = 0; i < nPlt; i++)
		{
			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
				Plants[i].type == "wind_offshore") {
				continue;
			}
			ex_fix += Plants[i].fix_cost * EV::Xop[n][i];
		}
		// var+fuel costs of plants
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
					Plants[i].type == "wind_offshore") {
					continue;
				}
				// var cost
				ex_var += time_weight[t] * Plants[i].var_cost * EV::prod[n][t][i];

				// fuel price to be updated later (dollar per thousand cubic feet=MMBTu)
				// NG fuel is handle in the NG network for case 2 and 3
				//if (Plants[i].type == "dfo")
				//{
				//	ex_thermal_fuel += time_weight[t] * dfo_pric * Plants[i].heat_rate * EV::prod[n][t][i];

				//}
				//else if (Plants[i].type == "coal")
				//{
				//	ex_thermal_fuel += time_weight[t] * coal_price * Plants[i].heat_rate * EV::prod[n][t][i];
				//}
				if (Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
				{
					ex_thermal_fuel += time_weight[t] * nuclear_price * Plants[i].heat_rate * EV::prod[n][t][i];
				}
				else if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
				{
					if (Setting::Case == 1)
					{
						ex_thermal_fuel += time_weight[t] * NG_price * Plants[i].heat_rate * EV::prod[n][t][i];
					}
				}

				if (Plants[i].type == "ng" || Plants[i].type == "CT" ||
					Plants[i].type == "CC" || Plants[i].type == "CC-CCS" ||
					Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
				{
					// startup cost
					ex_startup += time_weight[t] * Plants[i].startup_cost * EV::Xup[n][t][i];
				}

			}

			// load curtailment cost
			ex_shedd += time_weight[t] * E_curt_cost * EV::curtE[n][t];
		}


		// storage cost
		for (int r = 0; r < neSt; r++)
		{
			double s1 = (double)std::pow(1.0 / (1 + WACC), battery_lifetime);
			double capco = WACC / (1 - s1);
			ex_elec_str += capco * Estorage[r].power_cost * EV::YeCD[n][r] + capco * Estorage[r].energy_cost * EV::YeLev[n][r];
			ex_elec_str += Estorage[r].pFOM * EV::YeCD[n][r] + Estorage[r].eFOM * EV::YeLev[n][r]; // fixed cost per kw per year
		}
	}

	// investment cost of transmission lines
	for (int b = 0; b < Branches.size(); b++)
	{
		double s1 = (double)std::pow(1.0 / (1 + WACC), trans_line_lifespan);
		double capco = WACC / (1 - s1);
		int fb = Branches[b].from_bus;
		int tb = Branches[b].to_bus;
		// find the index of tb in the adj_buses of the fb
		int tbi = std::find(Enodes[fb].adj_buses.begin(), Enodes[fb].adj_buses.end(), tb) - Enodes[fb].adj_buses.begin();
		ex_trans += capco * trans_unit_cost * Branches[b].maxFlow * Branches[b].length * EV::Ze[b];
	}
	exp_Eobj = ex_est + ex_decom + ex_fix + ex_startup + ex_emis + ex_var + ex_thermal_fuel + ex_shedd + ex_trans + ex_elec_str;

	Model.addConstr(EV::est_cost == ex_est);
	Model.addConstr(EV::est_trans_cost == ex_trans);
	Model.addConstr(EV::decom_cost == ex_decom);
	Model.addConstr(EV::fixed_cost == ex_fix);
	Model.addConstr(EV::var_cost == ex_var);
	Model.addConstr(EV::thermal_fuel_cost == ex_thermal_fuel);
	Model.addConstr(EV::dfo_coal_emis_cost == ex_emis);
	Model.addConstr(EV::shedding_cost == ex_shedd);
	Model.addConstr(EV::elec_storage_cost == ex_elec_str);
	Model.addConstr(EV::startup_cost == ex_startup);

	Model.addConstr(exp_Eobj == EV::e_system_cost);
	Model.setObjective(exp_Eobj + exp_NGobj, GRB_MINIMIZE);
#pragma endregion


#pragma region Electricity Network Constraints
	// C1: number of generation units at each node	
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			// is plant type i part of initial plants types of node n:
			bool isInd = std::find(Enodes[n].Init_plt_ind.begin(), Enodes[n].Init_plt_ind.end(), i) != Enodes[n].Init_plt_ind.end();
			if (isInd)
			{
				// find the index of plant in the set of Init_plants of the node
				string tp = pltType2sym[i];
				int ind1 = std::find(Enodes[n].Init_plt_types.begin(), Enodes[n].Init_plt_types.end(), tp) - Enodes[n].Init_plt_types.begin();
				Model.addConstr(EV::Xop[n][i] == Enodes[n].Init_plt_count[ind1] - EV::Xdec[n][i] + EV::Xest[n][i]);
			}
			else
			{
				Model.addConstr(EV::Xop[n][i] == -EV::Xdec[n][i] + EV::Xest[n][i]);
			}

			for (int t = 0; t < Te.size(); t++)
			{
				// (w/o UC version) keep this to avoid power generation from hydro-new, dfor, etc.
				Model.addConstr(EV::prod[n][t][i] <= Plants[i].Pmax * EV::Xop[n][i]);
				// (w/o UC version)
				//			if (t > 0)
	//			{
	//				//Model.addConstr(Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i] >= -EV::prod[n][t][i] + EV::prod[n][t - 1][i]);
	//				Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
	//				Model.addConstr(-EV::prod[n][t][i] + EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
	//			}
			}
		}
	}

	//C2, C3, C4, C5: UC,  production limit, ramping for thermal units (ng, CT, CC, CC-CCS, nuclear)	
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			for (int t = 0; t < Te.size(); t++)
			{
				if (Plants[i].type == "ng" || Plants[i].type == "CT" ||
					Plants[i].type == "CC" || Plants[i].type == "CC-CCS" ||
					Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
				{
					// UC 1
					if (t > 0)
					{
						Model.addConstr(EV::X[n][t][i] - EV::X[n][t - 1][i] == EV::Xup[n][t][i] - EV::Xdown[n][t][i]);
					}
					// UC 2
					Model.addConstr(EV::X[n][t][i] <= EV::Xop[n][i]);
					// output limit
					Model.addConstr(EV::prod[n][t][i] >= Plants[i].Pmin * Plants[i].Pmax * EV::X[n][t][i]);
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].Pmax * EV::X[n][t][i]);
					//ramping
					if (t > 0)
					{
						double mx1 = std::max(Plants[i].Pmin, Plants[i].rampU);
						Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * (EV::X[n][t][i] - EV::Xup[n][t][i]) + mx1 * Plants[i].Pmax * EV::Xup[n][t][i]);
						Model.addConstr(-EV::prod[n][t][i] + EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * (EV::X[n][t][i] - EV::Xup[n][t][i]) + mx1 * Plants[i].Pmax * EV::Xup[n][t][i]);
					}
				}
			}
		}
	}

	// C5, C6: flow limit for electricity
	for (int br = 0; br < nBr; br++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (Branches[br].is_exist == 1)
			{
				Model.addConstr(EV::flowE[br][t] <= Branches[br].maxFlow);
				Model.addConstr(-EV::flowE[br][t] <= Branches[br].maxFlow);
				//Model.addConstr(-Branches[br].maxFlow <= EV::flowE[br][t]);
			}
			else
			{
				Model.addConstr(EV::flowE[br][t] <= Branches[br].maxFlow * EV::Ze[br]);
				Model.addConstr(-EV::flowE[br][t] <= Branches[br].maxFlow * EV::Ze[br]);
				//Model.addConstr(-Branches[br].maxFlow * EV::Ze[br] <= EV::flowE[br][t]);
			}
		}
	}

	// peak demand
	double max_dem = 0;
	for (int t = 0; t < Te.size(); t++)
	{
		double demm = 0;
		for (int n = 0; n < nEnode; n++)
		{
			demm += Enodes[n].demand[Te[t]];
		}
		if (demm > max_dem)
		{
			max_dem = demm;
		}
	}
	//C7: power balance
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			GRBLinExpr exp_prod(0);
			for (int i = 0; i < nPlt; i++)
			{
				exp_prod += EV::prod[n][t][i];
			}
			GRBLinExpr exp_trans(0);
			for (int m : Enodes[n].adj_buses)
			{
				// instead of defining a dictionray, one could traverse the entire branches...
				int key1 = n * 200 + m;
				// check if this key exist, if not, the other order exisit
				if (Le.count(key1) == 0)
				{
					key1 = m * 200 + n;
				}
				//each Le can contain multiple lines
				for (int l2 : Le[key1])
				{
					if (n > m)
					{
						exp_trans -= EV::flowE[l2][t];
					}
					else
					{
						exp_trans += EV::flowE[l2][t];
					}

				}
			}
			GRBLinExpr ex_store(0);
			for (int r = 0; r < neSt; r++)
			{
				ex_store += EV::eSdis[n][t][r] - EV::eSch[n][t][r];
			}
			double dem = Enodes[n].demand[Te[t]];


			// ignore trans
			//Model.addConstr(exp_prod + ex_store + EV::curtE[n][t] == dem);
			//Model.addConstr(exp_prod + exp_trans +  EV::curtE[n][t] == dem);

			Model.addConstr(exp_prod + exp_trans + ex_store + EV::curtE[n][t] == dem); //uncomment this line
		}
	}
	// C8: flow equation
	int ebr = 0;
	for (int br = 0; br < Branches.size(); br++)
	{
		int fb = Branches[br].from_bus;
		int tb = Branches[br].to_bus;
		for (int t = 0; t < Te.size(); t++)
		{
			//if (Branches[br].is_exist == 1 && !Setting::heuristics1_active)
			if (Branches[br].is_exist == 1)
			{
				Model.addConstr(EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) == 0);
				//Model.addConstr(-EV::flowE[br][t] + Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) <= 1e-4);
				//Model.addConstr(-1e-2 <= EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]));
			}
			//else if (!Setting::heuristics1_active)
			else
			{
				double Big_M = std::abs(Branches[br].suscep * 24);
				Model.addConstr(EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) <= Big_M * (1 - EV::Ze[br]));
				Model.addConstr(-EV::flowE[br][t] + Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) <= Big_M * (1 - EV::Ze[br]));
				//Model.addConstr(-Big_M * (1 - EV::Ze[br]) <= EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]));
			}
		}
	}

	// C9: phase angle (theta) limits. already applied in the definition of the variable
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 1; t < Te.size(); t++)
		{
			Model.addConstr(EV::theta[n][t] <= pi);
			Model.addConstr(-EV::theta[n][t] <= pi);
			//Model.addConstr(-pi <= EV::theta[n][t]);
		}
	}

	// C9.5: phase angle for reference node (node 0) is zero
	for (int t = 0; t < Te.size(); t++)
	{
		Model.addConstr(EV::theta[0][t] == 0);
	}

	// C10: VRE production profile
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "solar" || Plants[i].type == "solar-UPV")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::Xop[n][i]);
				}
				if (Plants[i].type == "wind" || Plants[i].type == "wind-new")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::Xop[n][i]);
				}
				if (Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new")
				{
					Model.addConstr(EV::prod[n][t][i] <= Plants[i].zonal_profile[Te[t]][0] * Plants[i].Pmax * EV::Xop[n][i]);
				}
				/*if (Plants[i].type == "hydro"|| Plants[i].type == "hydro-new")
				{
					Model.addConstr(EV::prod[n][t][i] <=Plants[i].Pmax * EV::Xop[n][i]);
				}*/
			}
		}
	}


	// C11: demand curtainlment constraint
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			double dem = Enodes[n].demand[Te[t]];
			//dem = dem/1.896;
			Model.addConstr(EV::curtE[n][t] <= dem);
		}
	}

	//C12: RPS constraints
	exp3 = 0;
	total_demand = 0;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			total_demand += time_weight[t] * Enodes[n].demand[Te[t]];

			for (int i = 0; i < nPlt; i++)
			{
				//exp2 += time_weight[t] * Plants[i].emis_rate * EV::prod[n][t][i];
				if (Plants[i].type == "solar" || Plants[i].type == "wind" ||
					Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new" ||
					Plants[i].type == "hydro" || Plants[i].type == "solar-UPV" ||
					Plants[i].type == "wind-new" || Plants[i].type == "hydro-new")
				{
					exp3 += time_weight[t] * EV::prod[n][t][i]; // production from VRE
				}
			}
		}
	}
	Model.addConstr(exp3 >= Setting::RPS * total_demand);

	// C14,C15,C16 storage constraints
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int r = 0; r < neSt; r++)
			{
				/*if (t == Te.size() - 1)
				{
					Model.addConstr(EV::eSlev[n][t][r] <= EV::eSlev[n][0][r]);
				}*/
				if (t == 0)
				{
					Model.addConstr(EV::eSlev[n][t][r] ==
						Estorage[r].eff_ch * EV::eSch[n][t][r] -
						(EV::eSdis[n][t][r] / Estorage[r].eff_disCh));
				}
				else
				{
					Model.addConstr(EV::eSlev[n][t][r] == EV::eSlev[n][t - 1][r] +
						Estorage[r].eff_ch * EV::eSch[n][t][r] -
						(EV::eSdis[n][t][r] / Estorage[r].eff_disCh));
				}
				Model.addConstr(EV::YeCD[n][r] >= EV::eSdis[n][t][r]);
				Model.addConstr(EV::YeCD[n][r] >= EV::eSch[n][t][r]);
				Model.addConstr(EV::YeLev[n][r] >= EV::eSlev[n][t][r]);
			}
		}
	}

	// start and ending of storage should be the same in case of using rep. days
	for (int n = 0; n < nEnode; n++)
	{
		for (int r = 0; r < neSt; r++)
		{
			for (int k = 0; k < Tg.size(); k++)
			{
				int str = k * 24;
				int end = (k + 1) * 24 - 1;
				Model.addConstr(EV::eSlev[n][str][r] == EV::eSlev[n][end][r]);
			}
		}
	}


	// C17: if an storage established or not
	for (int n = 0; n < nEnode; n++)
	{
		for (int r = 0; r < neSt; r++)
		{
			Model.addConstr(EV::YeLev[n][r] <= 10e8 * EV::YeStr[n][r]);
			Model.addConstr(EV::YeCD[n][r] <= EV::YeLev[n][r]);
		}
	}



#pragma endregion


#pragma region Solve and get values for UC vars
	Model.set(GRB_IntParam_OutputFlag, 0);
	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap * 3);
	Model.optimize();
	std::cout << "\n Second phase obj: " << Model.get(GRB_DoubleAttr_ObjVal);
	EV::val_X = new double** [nEnode];
	EV::val_Xup = new double** [nEnode];
	EV::val_Xdown = new double** [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_X[n] = new double* [Te.size()];
		EV::val_Xup[n] = new double* [Te.size()];
		EV::val_Xdown[n] = new double* [Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_X[n][t] = new double[nPlt]();
			EV::val_Xup[n][t] = new double[nPlt]();
			EV::val_Xdown[n][t] = new double[nPlt]();
			for (int i = 0; i < nPlt; i++)
			{
				EV::val_X[n][t][i] = EV::X[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_Xup[n][t][i] = EV::Xup[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_Xdown[n][t][i] = EV::Xdown[n][t][i].get(GRB_DoubleAttr_X);
			}
		}
	}
#pragma endregion





#pragma endregion

}



void SP_flow_Upper()
{
	// a heuristic to get the upper bound for flow variables in 
	// the power balance equation, to be used to add a cut in case all operational
	// variables are in the SP

	double eps = 0.00;
	//auto start = chrono::high_resolution_clock::now();
#pragma region Fetch Data

	std::map<string, int> sym2pltType = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7},
		{"Ct",8},{"CC",9},{"CC-CCS",10},{"solar-UPV",11},{"wind-new",12},
		{"wind-offshore-new",13},{"hydro-new",14},{"nuclear-new",15} };
	std::map<int, string> pltType2sym = { {0,"ng"},{1,"dfo"},
{2,"solar"},{3,"wind"},{4,"wind_offshore"},{5,"hydro"},{6,"coal"},{7,"nuclear"} };
	vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Params::Enodes.size();
	int nPlt = (int)Params::Plants.size();
	int nBr = (int)Params::Branches.size();
	int neSt = (int)Params::Estorage.size();
	vector<int> Tg = Params::Tg;
	vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	double pi = 3.141592;
	int nGnode = (int)Params::Gnodes.size();
	double WACC = Params::WACC;
	int trans_unit_cost = Params::trans_unit_cost;
	int trans_line_lifespan = Params::trans_line_lifespan;
	double NG_price = Params::NG_price;
	double dfo_pric = Params::dfo_pric;
	double coal_price = Params::coal_price;
	double nuclear_price = Params::nuclear_price;
	double E_curt_cost = Params::E_curt_cost;
	double G_curt_cost = Params::G_curt_cost;
	double pipe_per_mile = Params::pipe_per_mile;
	int pipe_lifespan = Params::pipe_lifespan;
	int battery_lifetime = Params::battery_lifetime;
	map<int, vector<int>> Le = Params::Le;
	vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion
	auto start = chrono::high_resolution_clock::now();
	GRBEnv* envB = 0;
	envB = new GRBEnv();
	GRBModel Model = GRBModel(envB);
	//GRBLinExpr exp_NGobj(0);
	GRBLinExpr exp_Eobj(0);

	Setting::relax_int_vars = true;
	Populate_EV_SP(Model);
	//Populate_GV(Model);
	//Elec_Module_Primal_SP(Model, exp_Eobj);

	//NG_Module(Model, exp_NGobj);

	Setting::Approach_1_active = false; Setting::Approach_2_active = false;

#pragma region Set some variables
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			//Model.addConstr(EV::flowE[b][t] == EV::val_flowE[b][t]);
		}
	}
#pragma endregion

	GRBVar** flowE_abs1 = new GRBVar * [nBr];
	GRBVar** flowE_abs2 = new GRBVar * [nBr];
#pragma region Objective Function
	for (int b = 0; b < nBr; b++)
	{
		flowE_abs1[b] = Model.addVars(Te.size());
		flowE_abs2[b] = Model.addVars(Te.size());
		for (int t = 0; t < Te.size(); t++)
		{
			exp_Eobj += flowE_abs1[b][t];
			//Model.addConstr(EV::flowE[b][t] == flowE_abs1[b][t] - flowE_abs2[b][t]);
			Model.addGenConstrAbs(flowE_abs1[b][t], EV::flowE[b][t]);
		}
	}
	Model.setObjective(exp_Eobj, GRB_MAXIMIZE);
	Model.addConstr(exp_Eobj <= 10e11);
#pragma endregion

#pragma region Electricity Network Constraints
	// C1, C2: number of generation units at each node	
	// C 3.5: max yearly generation for existing plants

	//C3, C4: production limit, ramping	
	// C5, C6: flow limit for electricity
	for (int br = 0; br < nBr; br++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (Branches[br].is_exist == 1)
			{
				Model.addConstr(-EV::flowE[br][t] >= -Branches[br].maxFlow);
				Model.addConstr(EV::flowE[br][t] >= -Branches[br].maxFlow);
			}
			else
			{
				Model.addConstr(-EV::flowE[br][t] >= -(Branches[br].maxFlow * EV::Ze[br]));
				Model.addConstr(EV::flowE[br][t] >= -(Branches[br].maxFlow * EV::Ze[br]));
			}
		}
	}

	// peak demand
	//C7: power balance
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			GRBLinExpr exp_prod(0);
			for (int i = 0; i < nPlt; i++)
			{
				exp_prod += EV::prod[n][t][i];
			}
			GRBLinExpr exp_trans(0);
			//double trans = 0;
			for (int m : Params::Enodes[n].adj_buses)
			{
				// instead of defining a dictionray, one could traverse all branches...
				int key1 = n * 200 + m;
				// check if this key exist, if not, the other order exisit
				if (Params::Le.count(key1) == 0)
				{
					key1 = m * 200 + n;
				}
				//each Le can contain multiple lines
				for (int l2 : Params::Le[key1])
				{
					if (n > m)
					{
						exp_trans -= EV::flowE[l2][t];
						//trans -= EV::val_flowE[l2][t];
					}
					else
					{
						exp_trans += EV::flowE[l2][t];
						//trans += EV::val_flowE[l2][t];
					}
				}
			}
			GRBLinExpr str(0);
			for (int r = 0; r < neSt; r++)
			{
				str += EV::eSdis[n][t][r] - EV::eSch[n][t][r]; //commented
			}
			double dem = Params::Enodes[n].demand[Te[t]];
			double rhs = dem;

			Model.addConstr(exp_trans + exp_prod + str + EV::curtE[n][t] == rhs);
		}
	}
	// C8: flow equation
	int ebr = 0;
	for (int br = 0; br < Branches.size(); br++)
	{
		int fb = Branches[br].from_bus;
		int tb = Branches[br].to_bus;
		for (int t = 0; t < Te.size(); t++)
		{
			//if (Branches[br].is_exist == 1 && !Setting::heuristics1_active)
			if (Branches[br].is_exist == 1)
			{
				Model.addConstr(EV::flowE[br][t] >= Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) - eps);
				Model.addConstr(-EV::flowE[br][t] >= -Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) - eps);
			}
			else
			{
				double Big_M = std::abs(Branches[br].suscep * 24);
				//double rhs = std::max(0.0, Big_M * (1 - EV::val_Ze[br]));
				Model.addConstr(-EV::flowE[br][t] + Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) >= -Big_M * (1 - EV::Ze[br]));
				Model.addConstr(EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) >= -Big_M * (1 - EV::Ze[br]));
				//Model.addConstr(-Big_M * (1 - EV::val_Ze[br]) <= EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]));
			}
		}
	}

	// C9: phase angle (theta) limits. already applied in the definition of the variable
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 1; t < Te.size(); t++)
		{
			Model.addConstr(-EV::theta[n][t] >= -pi);
			Model.addConstr(EV::theta[n][t] >= -pi);
		}
	}
	// C9.5: phase angle for reference node (node 0) is zero
	for (int t = 0; t < Te.size(); t++)
	{
		Model.addConstr(EV::theta[0][t] == 0);
	}

	// C10: VRE production profile


	// C11: demand curtainlment constraint
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			double dem = Enodes[n].demand[Te[t]];
			//dem = dem/1.896;
			Model.addConstr(EV::curtE[n][t] <= dem);
		}
	}
	////C12: RPS constraints
	//// C14,C15,C16 storage constraints


#pragma endregion


	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	//Model.set(GRB_IntParam_Method, 0);


	//Model.set(GRB_IntParam_OutputFlag, 0);

	Model.optimize();
	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize Primal SP!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
	}


	double obj_val = Model.get(GRB_DoubleAttr_ObjVal);
	double gap = 0;
	if (!Setting::relax_int_vars)
	{
		double gap = Model.get(GRB_DoubleAttr_MIPGap);
	}

	//Get_EV_vals(Model);
	//Print_Results(0, 0);

	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "\n\n flowE upper Elapsed time: " << Elapsed << endl;
	std::cout << "\t flowE_upper Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;




#pragma region Print variables


	std::cout << "val_flowE" << endl;
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			double flowe = flowE_abs1[b][t].get(GRB_DoubleAttr_X);
			if (flowe > 0)
			{
				std::cout << "flowE[" << b << "][" << t << "] = " << flowe << endl;
			}
		}
	}
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				double prod = EV::prod[n][t][i].get(GRB_DoubleAttr_X);
				if (prod > 0)
				{
					std::cout << "prod[" << n << "][" << t << "][" << i << "] = " << prod << endl;
				}
			}
		}
	}

	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			double curtE = EV::curtE[n][t].get(GRB_DoubleAttr_X);
			if (curtE > 0)
			{
				std::cout << "curtE[" << n << "][" << t << "] = " << curtE << endl;
			}
		}
	}
#pragma endregion

	Setting::Approach_1_active = true; Setting::Approach_2_active = false;
}


// validate Benders SP by solving the primal Subproblem
int Primal_subproblem(vector<SP>& Cuts)
{
	SP sp_new;
	double eps = Setting::Num_rep_days / 10; double eps2 = 0.02;
	//auto start = chrono::high_resolution_clock::now();
#pragma region Fetch Data

	std::map<string, int> sym2pltType = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7},
		{"CT",8},{"CC",9},{"CC-CCS",10},{"solar-UPV",11},{"wind-new",12},
		{"wind-offshore-new",13},{"hydro-new",14},{"nuclear-new",15} };
	std::map<int, string> pltType2sym = { {0,"ng"},{1,"dfo"},
{2,"solar"},{3,"wind"},{4,"wind_offshore"},{5,"hydro"},{6,"coal"},{7,"nuclear"} };
	vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Params::Enodes.size();
	int nPlt = (int)Params::Plants.size();
	int nBr = (int)Params::Branches.size();
	int neSt = (int)Params::Estorage.size();
	vector<int> Tg = Params::Tg;
	vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	double pi = 3.141592;
	int nGnode = (int)Params::Gnodes.size();
	double WACC = Params::WACC;
	int trans_unit_cost = Params::trans_unit_cost;
	int trans_line_lifespan = Params::trans_line_lifespan;
	double NG_price = Params::NG_price;
	double dfo_pric = Params::dfo_pric;
	double coal_price = Params::coal_price;
	double nuclear_price = Params::nuclear_price;
	double E_curt_cost = Params::E_curt_cost;
	double G_curt_cost = Params::G_curt_cost;
	double pipe_per_mile = Params::pipe_per_mile;
	int pipe_lifespan = Params::pipe_lifespan;
	int battery_lifetime = Params::battery_lifetime;
	map<int, vector<int>> Le = Params::Le;
	vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion
	auto start = chrono::high_resolution_clock::now();
	GRBEnv* envB = 0;
	envB = new GRBEnv();
	GRBModel Model = GRBModel(envB);
	//GRBLinExpr exp_NGobj(0);
	GRBLinExpr exp_Eobj(0);

	Setting::relax_int_vars = true;
	Populate_EV_SP(Model);
	//Populate_GV(Model);
	//Elec_Module_Primal_SP(Model, exp_Eobj);

	//NG_Module(Model, exp_NGobj);

	Setting::Approach_1_active = false; Setting::Approach_2_active = false;

#pragma region Set some variables
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			//Model.addConstr(EV::flowE[b][t] == EV::val_flowE[b][t]);
		}
	}
#pragma endregion

#pragma region Objective Function
	//Model.addConstr(exp_Eobj == 11);
	exp_Eobj = 0;
	Model.setObjective(exp_Eobj, GRB_MINIMIZE);
#pragma endregion

#pragma region Electricity Network Constraints
	// C1, C2: number of generation units at each node	
	// C 1: max yearly generation for existing plants
	//C3, C4: production limit, ramping	
	// C5, C6: flow limit for electricity
	for (int br = 0; br < nBr; br++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (Branches[br].is_exist == 1)
			{
				SP::d_delta11[br][t] = Model.addConstr(-EV::flowE[br][t] >= -Branches[br].maxFlow - eps);
				SP::d_delta12[br][t] = Model.addConstr(EV::flowE[br][t] >= -Branches[br].maxFlow - eps);
			}
			else
			{
				double rhs = (Branches[br].maxFlow * EV::val_Ze[br]);
				SP::d_delta21[br][t] = Model.addConstr(-EV::flowE[br][t] >= -rhs - eps);
				SP::d_delta22[br][t] = Model.addConstr(EV::flowE[br][t] >= -rhs - eps);
			}
		}
	}

	// peak demand
	// 
	//C7: power balance
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			double exp_prod = 0;
			for (int i = 0; i < nPlt; i++)
			{
				exp_prod += EV::val_prod[n][t][i];
			}
			GRBLinExpr exp_trans(0);
			//double trans = 0;
			for (int m : Params::Enodes[n].adj_buses)
			{
				// instead of defining a dictionray, one could traverse all branches...
				int key1 = n * 200 + m;
				// check if this key exist, if not, the other order exisit
				if (Params::Le.count(key1) == 0)
				{
					key1 = m * 200 + n;
				}
				//each Le can contain multiple lines
				for (int l2 : Params::Le[key1])
				{
					if (n > m)
					{
						exp_trans -= EV::flowE[l2][t];
						//trans -= EV::val_flowE[l2][t];
					}
					else
					{
						exp_trans += EV::flowE[l2][t];
						//trans += EV::val_flowE[l2][t];
					}
				}
			}
			double str = 0;
			for (int r = 0; r < neSt; r++)
			{
				str += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r]; //commented
			}
			double dem = Params::Enodes[n].demand[Te[t]];
			double rhs = dem - str - exp_prod - EV::val_curtE[n][t];

			//double rhs = dem - str - exp_prod;

			SP::d_theta1[n][t] = Model.addConstr(exp_trans >= rhs - eps);
			SP::d_theta2[n][t] = Model.addConstr(-exp_trans >= -rhs - eps);
		}
	}
	// C8: flow equation
	int ebr = 0;
	for (int br = 0; br < Branches.size(); br++)
	{
		int fb = Branches[br].from_bus;
		int tb = Branches[br].to_bus;
		for (int t = 0; t < Te.size(); t++)
		{
			//if (Branches[br].is_exist == 1 && !Setting::heuristics1_active)
			if (Branches[br].is_exist == 1)
			{
				SP::d_zeta21[br][t] = Model.addConstr(EV::flowE[br][t] >= Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) - eps);
				SP::d_zeta22[br][t] = Model.addConstr(-EV::flowE[br][t] >= -Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) - eps);
			}
			else
			{
				double Big_M = std::abs(Branches[br].suscep * 24);
				double rhs = std::max(0.0, Big_M * (1 - EV::val_Ze[br]));
				SP::d_zeta11[br][t] = Model.addConstr(-EV::flowE[br][t] + Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) >= -rhs - eps);
				SP::d_zeta12[br][t] = Model.addConstr(EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) >= -rhs - eps);
				//Model.addConstr(-Big_M * (1 - EV::val_Ze[br]) <= EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]));
			}
		}
	}

	// C9: phase angle (theta) limits. already applied in the definition of the variable
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 1; t < Te.size(); t++)
		{
			SP::d_eta1[n][t] = Model.addConstr(-EV::theta[n][t] >= -pi);
			SP::d_eta2[n][t] = Model.addConstr(EV::theta[n][t] >= -pi);
		}
	}
	// C9.5: phase angle for reference node (node 0) is zero
	for (int t = 0; t < Te.size(); t++)
	{
		Model.addConstr(EV::theta[0][t] == 0);
	}

	// C10: VRE production profile


	//// C11: demand curtainlment constraint

	////C12: RPS constraints
	//// C14,C15,C16 storage constraints

	// force flow>eps
	/*for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			Model.addConstr(EV::flowE[b][t] >= eps2);
			Model.addConstr(EV::flowE[b][t] <= -eps2);
		}
	}*/
#pragma endregion


	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model.set(GRB_IntParam_InfUnbdInfo, 1);// if additional info for infeasible/unbounded models is desired
	//Model.set(GRB_IntParam_OutputFlag, 0);

	Model.set(GRB_IntParam_Presolve, 0); // turn off presolve to be able to get extreme rays/farkas duals
	Model.set(GRB_IntParam_Method, 1); // 1: dual simplex, 0:simplex
	Model.optimize();
	return Model.get(GRB_IntAttr_Status);
	//if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	//{
	//	sp_new.dual_val_delta11 = new double* [nBr];
	//	sp_new.dual_val_delta12 = new double* [nBr];
	//	sp_new.dual_val_delta21 = new double* [nBr];
	//	sp_new.dual_val_delta22 = new double* [nBr];
	//	sp_new.dual_val_theta1 = new double* [nEnode];
	//	sp_new.dual_val_theta2 = new double* [nEnode];
	//	sp_new.dual_val_zeta11 = new double* [nBr];
	//	sp_new.dual_val_zeta12 = new double* [nBr];
	//	sp_new.dual_val_zeta21 = new double* [nBr];
	//	sp_new.dual_val_zeta22 = new double* [nBr];
	//	sp_new.dual_val_eta1 = new double* [nEnode];
	//	sp_new.dual_val_eta2 = new double* [nEnode];
	//	double fd = 0;
	//	for (int n = 0; n < nEnode; n++)
	//	{
	//		sp_new.dual_val_theta1[n] = new double[Te.size()]();
	//		sp_new.dual_val_theta2[n] = new double[Te.size()]();

	//		sp_new.dual_val_eta1[n] = new double[Te.size()]();
	//		sp_new.dual_val_eta2[n] = new double[Te.size()]();
	//		double ave_price = 0;
	//		for (int t = 0; t < Te.size(); t++)
	//		{
	//			sp_new.dual_val_theta1[n][t] = std::max(0.0, sp_new.d_theta1[n][t].get(GRB_DoubleAttr_FarkasDual));
	//			sp_new.dual_val_theta2[n][t] = std::max(0.0, sp_new.d_theta2[n][t].get(GRB_DoubleAttr_FarkasDual));
	//			fd += std::abs(sp_new.dual_val_theta1[n][t]) + std::abs(sp_new.dual_val_theta2[n][t]);
	//			if (t > 0)
	//			{
	//				sp_new.dual_val_eta1[n][t] = std::max(0.0, sp_new.d_eta1[n][t].get(GRB_DoubleAttr_FarkasDual));
	//				fd += std::abs(sp_new.dual_val_eta1[n][t]);
	//				sp_new.dual_val_eta2[n][t] = std::max(0.0, sp_new.d_eta2[n][t].get(GRB_DoubleAttr_FarkasDual));
	//				fd += std::abs(sp_new.dual_val_eta2[n][t]);
	//			}
	//		}
	//	}

	//	for (int b = 0; b < nBr; b++)
	//	{
	//		sp_new.dual_val_delta11[b] = new double[Te.size()]();
	//		sp_new.dual_val_delta12[b] = new double[Te.size()]();
	//		sp_new.dual_val_delta21[b] = new double[Te.size()]();
	//		sp_new.dual_val_delta22[b] = new double[Te.size()]();
	//		sp_new.dual_val_zeta11[b] = new double[Te.size()]();
	//		sp_new.dual_val_zeta12[b] = new double[Te.size()]();
	//		sp_new.dual_val_zeta21[b] = new double[Te.size()]();
	//		sp_new.dual_val_zeta22[b] = new double[Te.size()]();
	//		for (int t = 0; t < Te.size(); t++)
	//		{
	//			if (Branches[b].is_exist == 1)
	//			{
	//				sp_new.dual_val_delta11[b][t] = std::max(0.0, SP::d_delta11[b][t].get(GRB_DoubleAttr_FarkasDual));
	//				sp_new.dual_val_delta12[b][t] = std::max(0.0, SP::d_delta12[b][t].get(GRB_DoubleAttr_FarkasDual));
	//				sp_new.dual_val_zeta21[b][t] = std::max(0.0, SP::d_zeta21[b][t].get(GRB_DoubleAttr_FarkasDual));
	//				sp_new.dual_val_zeta22[b][t] = std::max(0.0, SP::d_zeta22[b][t].get(GRB_DoubleAttr_FarkasDual));
	//			}
	//			else
	//			{
	//				sp_new.dual_val_delta21[b][t] = std::max(0.0, SP::d_delta21[b][t].get(GRB_DoubleAttr_FarkasDual));
	//				sp_new.dual_val_delta22[b][t] = std::max(0.0, SP::d_delta22[b][t].get(GRB_DoubleAttr_FarkasDual));


	//				sp_new.dual_val_zeta11[b][t] = std::max(0.0, SP::d_zeta11[b][t].get(GRB_DoubleAttr_FarkasDual));
	//				sp_new.dual_val_zeta12[b][t] = std::max(0.0, SP::d_zeta12[b][t].get(GRB_DoubleAttr_FarkasDual));
	//			}
	//			//fd += SP::dual_val_delta11[b][t] + SP::dual_val_delta12[b][t] + SP::dual_val_delta21[b][t] + SP::dual_val_delta22[b][t];

	//		}
	//	}

	//	double fd2 = SP::d_theta1[0][2].get(GRB_DoubleAttr_FarkasDual);


	//	//Model.get(GRB_DoubleAttr_FarkasDual,2, 3);
	//	std::cout << "Failed to optimize Primal SP!!!" << endl;
	//	std::cout << Model.get(GRB_IntAttr_Status);
	//}


	//double obj_val = Model.get(GRB_DoubleAttr_ObjVal);
	//double gap = 0;
	//if (!Setting::relax_int_vars)
	//{
	//	double gap = Model.get(GRB_DoubleAttr_MIPGap);
	//}

	////Get_EV_vals(Model);
	////Print_Results(0, 0);

	//auto end = chrono::high_resolution_clock::now();
	//double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	//std::cout << "\n\n Primal SP Elapsed time: " << Elapsed << endl;
	//std::cout << "\t Primal SP Obj Value:" << obj_val << endl;
	//std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;


#pragma region SP primal objective value in IM
	//GRBLinExpr ex_var(0);
	//GRBLinExpr ex_thermal_fuel(0);

	//GRBLinExpr ex_shedd(0);
	//SP::SP_Primal_obj = 0;
	//for (int n = 0; n < nEnode; n++)
	//{

	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		for (int i = 0; i < nPlt; i++)
	//		{
	//			if (Plants[i].type == "dfo" || Plants[i].type == "coal" ||
	//				Plants[i].type == "wind_offshore") {
	//				continue;
	//			}
	//			// var cost
	//			SP::SP_Primal_obj += time_weight[t] * Plants[i].var_cost * EV::val_prod[n][t][i];

	//			// fuel price to be updated later (dollar per thousand cubic feet=MMBTu)
	//			// NG fuel is handle in the NG network for case 2 and 3
	//			if (Plants[i].type == "dfo")
	//			{
	//				SP::SP_Primal_obj += time_weight[t] * dfo_pric * Plants[i].heat_rate * EV::val_prod[n][t][i];

	//			}
	//			else if (Plants[i].type == "coal")
	//			{
	//				SP::SP_Primal_obj += time_weight[t] * coal_price * Plants[i].heat_rate * EV::val_prod[n][t][i];
	//			}
	//			else if (Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
	//			{
	//				SP::SP_Primal_obj += time_weight[t] * nuclear_price * Plants[i].heat_rate * EV::val_prod[n][t][i];
	//			}
	//			else if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
	//			{
	//				if (Setting::Case == 1)
	//				{
	//					SP::SP_Primal_obj += time_weight[t] * NG_price * Plants[i].heat_rate * EV::val_prod[n][t][i];
	//				}
	//			}

	//			// emission cost (no longer needed)
	//			/*if (Plants[i].type != "ng" && Plants[i].type != "CT" && Plants[i].type != "CC" && Plants[i].type != "CC-CCS")
	//			{
	//				ex_emis += time_weight[t] * Plants[i].emis_cost * Plants[i].emis_rate * EV::prod[n][t][i];
	//			}*/
	//		}

	//		// load curtailment cost
	//		SP::SP_Primal_obj += time_weight[t] * E_curt_cost * EV::val_curtE[n][t];
	//	}


	//	// storage cost

	//}
	//std::cout << "\n\n \t\t Primal problem cost in IM :" << SP::SP_Primal_obj << endl;
#pragma endregion

#pragma region Get dual values
	/*sp_new.dual_val_delta11 = new double* [nBr];
	sp_new.dual_val_delta12 = new double* [nBr];
	sp_new.dual_val_delta21 = new double* [nBr];
	sp_new.dual_val_delta22 = new double* [nBr];
	sp_new.dual_val_theta1 = new double* [nEnode];
	sp_new.dual_val_theta2 = new double* [nEnode];
	sp_new.dual_val_zeta11 = new double* [nBr];
	sp_new.dual_val_zeta12 = new double* [nBr];
	sp_new.dual_val_zeta21 = new double* [nBr];
	sp_new.dual_val_zeta22 = new double* [nBr];
	sp_new.dual_val_eta1 = new double* [nEnode];
	sp_new.dual_val_eta2 = new double* [nEnode];

	for (int n = 0; n < nEnode; n++)
	{
		sp_new.dual_val_theta1[n] = new double[Te.size()]();
		sp_new.dual_val_theta2[n] = new double[Te.size()]();

		sp_new.dual_val_eta1[n] = new double[Te.size()]();
		sp_new.dual_val_eta2[n] = new double[Te.size()]();
		double ave_price = 0;
		for (int t = 0; t < Te.size(); t++)
		{
			sp_new.dual_val_theta1[n][t] = std::max(0.0, sp_new.d_theta1[n][t].get(GRB_DoubleAttr_Pi));
			sp_new.dual_val_theta2[n][t] = std::max(0.0, sp_new.d_theta2[n][t].get(GRB_DoubleAttr_Pi));
			if (t > 0)
			{
				sp_new.dual_val_eta1[n][t] = std::max(0.0, sp_new.d_eta1[n][t].get(GRB_DoubleAttr_Pi));
				sp_new.dual_val_eta2[n][t] = std::max(0.0, sp_new.d_eta2[n][t].get(GRB_DoubleAttr_Pi));
			}
		}
	}

	for (int b = 0; b < nBr; b++)
	{
		sp_new.dual_val_delta11[b] = new double[Te.size()]();
		sp_new.dual_val_delta12[b] = new double[Te.size()]();
		sp_new.dual_val_delta21[b] = new double[Te.size()]();
		sp_new.dual_val_delta22[b] = new double[Te.size()]();
		sp_new.dual_val_zeta11[b] = new double[Te.size()]();
		sp_new.dual_val_zeta12[b] = new double[Te.size()]();
		sp_new.dual_val_zeta21[b] = new double[Te.size()]();
		sp_new.dual_val_zeta22[b] = new double[Te.size()]();
		for (int t = 0; t < Te.size(); t++)
		{
			if (Branches[b].is_exist == 1)
			{
				sp_new.dual_val_delta11[b][t] = std::max(0.0, SP::d_delta11[b][t].get(GRB_DoubleAttr_Pi));
				sp_new.dual_val_delta12[b][t] = std::max(0.0, SP::d_delta12[b][t].get(GRB_DoubleAttr_Pi));
				sp_new.dual_val_zeta21[b][t] = std::max(0.0, SP::d_zeta21[b][t].get(GRB_DoubleAttr_Pi));
				sp_new.dual_val_zeta22[b][t] = std::max(0.0, SP::d_zeta22[b][t].get(GRB_DoubleAttr_Pi));
			}
			else
			{
				sp_new.dual_val_delta21[b][t] = std::max(0.0, SP::d_delta21[b][t].get(GRB_DoubleAttr_Pi));
				sp_new.dual_val_delta22[b][t] = std::max(0.0, SP::d_delta22[b][t].get(GRB_DoubleAttr_Pi));


				sp_new.dual_val_zeta11[b][t] = std::max(0.0, SP::d_zeta11[b][t].get(GRB_DoubleAttr_Pi));
				sp_new.dual_val_zeta12[b][t] = std::max(0.0, SP::d_zeta12[b][t].get(GRB_DoubleAttr_Pi));
			}
		}
	}

	std::cout << "dual_val_theta1: " << endl;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (sp_new.dual_val_theta1[n][t] > 0.0)
			{
				std::cout << sp_new.dual_val_theta1[n][t] << endl;
			}
		}
	}
	std::cout << "dual_val_theta2: " << endl;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (sp_new.dual_val_theta2[n][t] > 0.0)
			{
				std::cout << sp_new.dual_val_theta2[n][t] << endl;
			}
		}
	}*/


	// coupling constraints 
	/*SP::dual_val_rho1 = new double** [nGnode];
	SP::dual_val_rho2 = new double** [nGnode];
	SP::dual_val_tau = SP::d_tau.get(GRB_DoubleAttr_Pi);
	for (int k = 0; k < nGnode; k++)
	{
		SP::dual_val_rho1[k] = new double* [nEnode];
		SP::dual_val_rho2[k] = new double* [nEnode];
		for (int n : Gnodes[k].adjE)
		{
			SP::dual_val_rho1[k][n] = new double[Tg.size()]();
			SP::dual_val_rho2[k][n] = new double[Tg.size()]();
			for (int tau = 0; tau < Tg.size(); tau++)
			{
				SP::dual_val_rho1[k][n][tau] = SP::d_rho1[k][n][tau].get(GRB_DoubleAttr_Pi);
				SP::dual_val_rho2[k][n][tau] = SP::d_rho2[k][n][tau].get(GRB_DoubleAttr_Pi);
			}
		}
	}*/
#pragma endregion

#pragma region calculate Dual SP objective value
	//SP::SP_Dual_obj = 0;
	//for (int i = 0; i < nPlt; i++)
	//{
	//	if (Plants[i].is_exis != 1) { continue; } // only applies to existing plants
	//	SP::SP_Dual_obj -= Plants[i].max_yearly_gen * SP::dual_val_alpha[i];
	//}

	//for (int n = 0; n < nEnode; n++)
	//{
	//	/*for (int i = 0; i < nPlt; i++)
	//	{
	//		if (EV::val_Xop[n][i] > 0)
	//		{
	//			std::cout << "xop[" << n << "][" << i << "] = " << EV::val_Xop[n][i] << endl;
	//		}
	//	}*/
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		for (int i = 0; i < nPlt; i++)
	//		{
	//			SP::SP_Dual_obj -= Plants[i].Pmax * EV::val_Xop[n][i] * SP::dual_val_beta[n][t][i];

	//			if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
	//			{
	//				if (t > 0)
	//				{
	//					SP::SP_Dual_obj -= Plants[i].rampU * Plants[i].Pmax * EV::val_Xop[n][i] * SP::dual_val_gamma1[n][t][i];
	//					SP::SP_Dual_obj -= Plants[i].rampU * Plants[i].Pmax * EV::val_Xop[n][i] * SP::dual_val_gamma2[n][t][i];
	//				}
	//			}

	//			if (Plants[i].type == "solar" || Plants[i].type == "solar-UPV")
	//			{
	//				SP::SP_Dual_obj -= (eps + Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::val_Xop[n][i]) * SP::dual_val_pi[n][t][i];
	//			}
	//			if (Plants[i].type == "wind" || Plants[i].type == "wind-new")
	//			{
	//				SP::SP_Dual_obj -= (eps + Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::val_Xop[n][i]) * SP::dual_val_pi[n][t][i];
	//			}
	//			if (Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new")
	//			{
	//				SP::SP_Dual_obj -= (eps + Plants[i].zonal_profile[Te[t]][0] * Plants[i].Pmax * EV::val_Xop[n][i]) * SP::dual_val_pi[n][t][i];
	//			}
	//		}
	//		double dem = Enodes[n].demand[Te[t]];
	//		double str = 0;
	//		for (int r = 0; r < neSt; r++)
	//		{
	//			str += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r];
	//		}
	//		SP::SP_Dual_obj += (dem - str) * SP::dual_val_theta1[n][t]; // commented
	//		SP::SP_Dual_obj -= (dem - str) * SP::dual_val_theta2[n][t];

	//		SP::SP_Dual_obj -= dem * SP::dual_val_phi[n][t];
	//		SP::SP_Dual_obj += Setting::RPS * time_weight[t] * dem * SP::dual_val_omega; //commented

	//		SP::SP_Dual_obj -= pi * SP::dual_val_eta1[n][t];
	//		SP::SP_Dual_obj -= pi * SP::dual_val_eta2[n][t];
	//	}
	//}

	//for (int b = 0; b < nBr; b++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (Branches[b].is_exist == 1)
	//		{
	//			SP::SP_Dual_obj -= Branches[b].maxFlow * SP::dual_val_delta11[b][t];
	//			SP::SP_Dual_obj -= Branches[b].maxFlow * SP::dual_val_delta12[b][t];
	//		}
	//		else
	//		{
	//			double Big_M = std::abs(Branches[b].suscep * 24);
	//			SP::SP_Dual_obj -= Branches[b].maxFlow * EV::val_Ze[b] * SP::dual_val_delta21[b][t];
	//			SP::SP_Dual_obj -= Branches[b].maxFlow * EV::val_Ze[b] * SP::dual_val_delta22[b][t];

	//			SP::SP_Dual_obj += Big_M * (1 - EV::val_Ze[b]) * SP::dual_val_zeta11[b][t];
	//			SP::SP_Dual_obj -= Big_M * (1 - EV::val_Ze[b]) * SP::dual_val_zeta12[b][t];
	//		}
	//	}
	//}

	//double rem = Setting::Emis_redu_goal * Setting::PE;

	//for (int k = 0; k < nGnode; k++)
	//{
	//	for (int tau = 0; tau < Tg.size(); tau++)
	//	{
	//		double fge = 0;
	//		for (int n : Gnodes[k].adjE)
	//		{
	//			SP::SP_Dual_obj += GV::val_flowGE[k][n][tau] * SP::dual_val_rho1[k][n][tau];
	//			SP::SP_Dual_obj -= GV::val_flowGE[k][n][tau] * SP::dual_val_rho2[k][n][tau];
	//			//std::cout << SP::dual_val_rho[k][n][tau] << endl;
	//			//std::cout << "flowGE[" << k << "][" << n << "][" << tau << "] = " << GV::val_flowGE[k][n][tau] << endl;
	//			fge += GV::val_flowGE[k][n][tau];
	//		}
	//		rem -= (RepDaysCount[tau] * Params::NG_emis_rate * (GV::val_supply[k][tau] - fge));//commented
	//	}
	//}
	//SP::SP_Dual_obj -= rem * SP::dual_val_tau;
	//std::cout << "\t\t Dual SP objective value: " << SP::SP_Dual_obj << endl;
#pragma endregion


#pragma region demo calc
	//double null_obj = 0;
	//for (int i = 0; i < nPlt; i++)
	//{
	//	if (Plants[i].is_exis != 1) { continue; } // only applies to existing plants
	//	null_obj += Plants[i].max_yearly_gen * SP::dual_val_alpha[i];
	//}

	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{

	//		double dem = Enodes[n].demand[Te[t]];
	//		double str = 0;
	//		for (int r = 0; r < neSt; r++)
	//		{
	//			str += 0;
	//		}
	//		null_obj += (dem - str) * SP::dual_val_theta1[n][t]; // commented

	//		null_obj += dem * SP::dual_val_phi[n][t];
	//		null_obj += Setting::RPS * time_weight[t] * dem * SP::dual_val_omega; //commented

	//		SP::SP_Dual_obj += pi * SP::dual_val_eta1[n][t] + pi * SP::dual_val_eta2[n][t];
	//	}
	//}

	//for (int b = 0; b < nBr; b++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (Branches[b].is_exist == 1)
	//		{
	//			null_obj += Branches[b].maxFlow * (SP::dual_val_delta11[b][t] + SP::dual_val_delta12[b][t]);
	//		}
	//		else
	//		{
	//			double Big_M = std::abs(Branches[b].suscep * 24);
	//			null_obj += Big_M * (1) * (SP::dual_val_zeta11[b][t] + SP::dual_val_zeta12[b][t]);
	//		}
	//	}
	//}

	//rem = Setting::Emis_redu_goal * Setting::PE;

	//for (int k = 0; k < nGnode; k++)
	//{
	//	for (int tau = 0; tau < Tg.size(); tau++)
	//	{
	//		double fge = 0;
	//		for (int n : Gnodes[k].adjE)
	//		{
	//			//SP::SP_Dual_obj += GV::val_flowGE[k][n][tau] * SP::dual_val_rho[k][n][tau]; // commented
	//			//std::cout << SP::dual_val_rho[k][n][tau] << endl;
	//			//std::cout << "flowGE[" << k << "][" << n << "][" << tau << "] = " << GV::val_flowGE[k][n][tau] << endl;
	//			//fge += GV::val_flowGE[k][n][tau];
	//		}
	//		rem -= (RepDaysCount[tau] * Params::NG_emis_rate * (GV::val_supply[k][tau] - 0));//commented
	//	}
	//}
	//null_obj += rem * SP::dual_val_tau;
	//std::cout << "\t\t Null objective value: " << null_obj << endl;
#pragma endregion



#pragma region Print variables
	//std::cout << "dual_val_theta1: " << endl;
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (sp_new.dual_val_theta1[n][t] > 0.0)
	//		{
	//			std::cout << sp_new.dual_val_theta1[n][t] << endl;
	//		}
	//	}
	//}
	//std::cout << "dual_val_theta2: " << endl;
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (sp_new.dual_val_theta2[n][t] > 0.0)
	//		{
	//			std::cout << sp_new.dual_val_theta2[n][t] << endl;
	//		}
	//	}
	//}

	//std::cout << "dual_val_eta1: " << endl;
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (sp_new.dual_val_eta1[n][t] > 0.0)
	//		{
	//			std::cout << sp_new.dual_val_eta1[n][t] << endl;
	//		}
	//	}
	//}
	//std::cout << "dual_val_eta2: " << endl;
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		if (sp_new.dual_val_eta2[n][t] > 0.0)
	//		{
	//			std::cout << sp_new.dual_val_eta2[n][t] << endl;
	//		}
	//	}
	//}
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		EV::val_theta[n][t] = EV::theta[n][t].get(GRB_DoubleAttr_X);
	//		if (EV::val_theta[n][t] > 0.0)
	//		{
	//			std::cout << EV::theta[n][t].get(GRB_DoubleAttr_X) << endl;
	//		}
	//	}
	//}

	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		double exp_prod = 0;
	//		for (int i = 0; i < nPlt; i++)
	//		{
	//			exp_prod += EV::val_prod[n][t][i];
	//		}
	//		double exp_trans = 0;
	//		//double trans = 0;
	//		for (int m : Params::Enodes[n].adj_buses)
	//		{
	//			// instead of defining a dictionray, one could traverse all branches...
	//			int key1 = n * 200 + m;
	//			// check if this key exist, if not, the other order exisit
	//			if (Params::Le.count(key1) == 0)
	//			{
	//				key1 = m * 200 + n;
	//			}
	//			//each Le can contain multiple lines
	//			for (int l2 : Params::Le[key1])
	//			{
	//				if (n > m)
	//				{
	//					exp_trans -= EV::flowE[l2][t].get(GRB_DoubleAttr_X);
	//					//trans -= EV::val_flowE[l2][t];
	//				}
	//				else
	//				{
	//					exp_trans += EV::flowE[l2][t].get(GRB_DoubleAttr_X);
	//					//trans += EV::val_flowE[l2][t];
	//				}
	//			}
	//		}
	//		double str = 0;
	//		for (int r = 0; r < neSt; r++)
	//		{
	//			str += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r]; //commented
	//		}
	//		double dem = Params::Enodes[n].demand[Te[t]];
	//		//double rhs = dem - str - exp_prod - EV::val_curtE[n][t];

	//		double rhs = dem - str - exp_prod;
	//		cout << "lhs: " << exp_trans << "\t rhs: " << rhs << endl;;
	//	}
	//}


	/*std::cout << "val_flowE" << endl;
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_flowE[b][t] = EV::flowE[b][t].get(GRB_DoubleAttr_X);
			if (std::abs(EV::val_flowE[b][t]) > 0.01)
			{
				std::cout << "flowE[" << b << "][" << t << "] = " << EV::val_flowE[b][t] << endl;
			}
		}
	}*/


	//for (int br = 0; br < 10; br++)
	//{
	//	int fb = Branches[br].from_bus;
	//	int tb = Branches[br].to_bus;
	//	for (int t = 0; t < 3; t++)
	//	{
	//		if (Branches[br].is_exist == 1)
	//		{
	//			//Model.addConstr(EV::flowE[br][t] >= Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) - eps);
	//			//Model.addConstr(-EV::flowE[br][t] >= -Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) - eps);
	//			//std::cout << "theta[from]: " << EV::val_theta[fb][t] << "\t theta[to]: " << EV::val_theta[tb][t] << endl;
	//			//std::cout << "flowE: " << EV::val_flowE[br][t] << "\t rhs: " << Branches[br].suscep * (EV::val_theta[tb][t] - EV::val_theta[fb][t]) << endl;

	//		}
	//	}
	//}
#pragma endregion

	//Cuts.push_back(sp_new);
	Setting::Approach_1_active = true; Setting::Approach_2_active = false;
}

