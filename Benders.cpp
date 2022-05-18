#include"Models_Funcs.h";

#pragma region SP variables/value
GRBVar*** SP::rho;
GRBVar SP::tau;
GRBVar* SP::alpha;
GRBVar*** SP::beta;
GRBVar*** SP::gamma1;
GRBVar*** SP::gamma2;
GRBVar** SP::delta11;
GRBVar** SP::delta12;
GRBVar** SP::delta21;
GRBVar** SP::delta22;
GRBVar** SP::theta;
GRBVar** SP::zeta11;
GRBVar** SP::zeta12;
GRBVar** SP::zeta21;
GRBVar** SP::zeta22;
GRBVar** SP::eta1;
GRBVar** SP::eta2;
GRBVar* SP::eta3;
GRBVar*** SP::pi;
GRBVar** SP::phi;
GRBVar SP::omega;
#pragma endregion

// validate Benders SP by solving the primal Subproblem


void Elec_Module_Primal_SP(GRBModel& Model, GRBLinExpr& exp_Eobj)
{
	//auto start = chrono::high_resolution_clock::now();
#pragma region Fetch Data

	std::map<string, int> sym2pltType = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7},
		{"Ct",8},{"CC",9},{"CC-CCS",10},{"solar-UPV",11},{"wind-new",12},
		{"wind-offshore-new",13},{"hydro-new",14},{"nuclear-new",15} };
	std::map<int, string> pltType2sym = { {0,"ng"},{1,"dfo"},
{2,"solar"},{3,"wind"},{4,"wind_offshore"},{5,"hydro"},{6,"coal"},{7,"nuclear"} };
	/*vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<branch> Branches = Params::Branches;*/
	int nEnode = (int)Params::Enodes.size();
	int nPlt = (int)Params::Plants.size();
	int nBr = (int)Params::Branches.size();
	int neSt = (int)Params::Estorage.size();
	//vector<int> Tg = Params::Tg;
	//vector<int> Te = Params::Te;
	//vector<int> time_weight = Params::time_weight;
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
	//map<int, vector<int>> Le = Params::Le;
	//vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion


#pragma region Objective Function
	GRBLinExpr ex_est(0);
	GRBLinExpr ex_decom(0);
	GRBLinExpr ex_fix(0);
	GRBLinExpr ex_var(0);
	GRBLinExpr ex_thermal_fuel(0);
	GRBLinExpr ex_emis(0);
	GRBLinExpr ex_shedd(0);
	GRBLinExpr ex_trans(0);
	GRBLinExpr ex_elec_str(0);


	for (int n = 0; n < nEnode; n++)
	{
		//for (int i = 0; i < nPlt; i++)
		//{
		//	// Investment/decommission cost of plants
		//	if (Plants[i].type == "dfo" || Plants[i].type == "coal" || Plants[i].type == "wind_offshore") { continue; }
		//	double s1 = (double)std::pow(1.0 / (1 + WACC), Plants[i].lifetime);
		//	double capco = WACC / (1 - s1) * Plants[i].Reg_coeffs_per_state[n];
		//	ex_est += capco * Plants[i].capex * EV::Xest[n][i];
		//	ex_decom += Plants[i].decom_cost * EV::Xdec[n][i];
		//}

		//// fixed cost (annual, so no iteration over time)
		//for (int i = 0; i < nPlt; i++)
		//{
		//	if (Plants[i].type == "dfo" || Plants[i].type == "coal" || Plants[i].type == "wind_offshore") { continue; }
		//	ex_fix += Plants[i].fix_cost * EV::Xop[n][i];
		//}
		// var+fuel costs of plants
		for (int t = 0; t < Params::Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Params::Plants[i].type == "dfo" || Params::Plants[i].type == "coal" ||
					Params::Plants[i].type == "wind_offshore") {
					continue;
				}
				// var cost
				ex_var += Params::time_weight[t] * Params::Plants[i].var_cost * EV::prod[n][t][i];

				// fuel price to be updated later (dollar per thousand cubic feet=MMBTu)
				// NG fuel is handle in the NG network for case 2 and 3
				if (Params::Plants[i].type == "dfo")
				{
					ex_thermal_fuel += Params::time_weight[t] * dfo_pric * Params::Plants[i].heat_rate * EV::prod[n][t][i];

				}
				else if (Params::Plants[i].type == "coal")
				{
					ex_thermal_fuel += Params::time_weight[t] * coal_price * Params::Plants[i].heat_rate * EV::prod[n][t][i];
				}
				else if (Params::Plants[i].type == "nuclear" || Params::Plants[i].type == "nuclear-new")
				{
					ex_thermal_fuel += Params::time_weight[t] * nuclear_price * Params::Plants[i].heat_rate * EV::prod[n][t][i];
				}
				else if (Params::Plants[i].type == "ng" || Params::Plants[i].type == "CT" || Params::Plants[i].type == "CC" || Params::Plants[i].type == "CC-CCS")
				{
					if (Setting::Case == 1)
					{
						ex_thermal_fuel += Params::time_weight[t] * NG_price * Params::Plants[i].heat_rate * EV::prod[n][t][i];
					}
				}

				// emission cost (no longer needed)
				/*if (Plants[i].type != "ng" && Plants[i].type != "CT" && Plants[i].type != "CC" && Plants[i].type != "CC-CCS")
				{
					ex_emis += time_weight[t] * Plants[i].emis_cost * Plants[i].emis_rate * EV::prod[n][t][i];
				}*/
			}

			// load curtailment cost
			ex_shedd += Params::time_weight[t] * E_curt_cost * EV::curtE[n][t];
		}


		// storage cost
		for (int r = 0; r < neSt; r++)
		{
			double s1 = (double)std::pow(1.0 / (1 + WACC), battery_lifetime);
			double capco = WACC / (1 - s1);
			ex_elec_str += capco * Params::Estorage[r].power_cost * EV::YeCD[n][r] + capco * Params::Estorage[r].energy_cost * EV::YeLev[n][r];
			ex_elec_str += Params::Estorage[r].pFOM * EV::YeCD[n][r] + Params::Estorage[r].eFOM * EV::YeLev[n][r]; // fixed cost per kw per year
		}
	}

	// investment cost of transmission lines
	//for (int b = 0; b < Branches.size(); b++)
	//{
	//	double s1 = (double)std::pow(1.0 / (1 + WACC), trans_line_lifespan);
	//	double capco = WACC / (1 - s1);
	//	int fb = Branches[b].from_bus;
	//	int tb = Branches[b].to_bus;
	//	// find the index of tb in the adj_buses of the fb
	//	int tbi = std::find(Enodes[fb].adj_buses.begin(), Enodes[fb].adj_buses.end(), tb) - Enodes[fb].adj_buses.begin();
	//	ex_trans += capco * trans_unit_cost * Branches[b].maxFlow * Branches[b].length * EV::Ze[b];
	//}
	exp_Eobj = ex_emis + ex_var + ex_thermal_fuel + ex_shedd + ex_elec_str;

	//Model.addConstr(EV::est_cost == ex_est);
	//Model.addConstr(EV::decom_cost == ex_decom);
	//Model.addConstr(EV::fixed_cost == ex_fix);
	Model.addConstr(EV::var_cost == ex_var);
	Model.addConstr(EV::thermal_fuel_cost == ex_thermal_fuel);
	Model.addConstr(EV::dfo_coal_emis_cost == ex_emis);
	Model.addConstr(EV::shedding_cost == ex_shedd);
	Model.addConstr(EV::elec_storage_cost == ex_elec_str);

	Model.addConstr(exp_Eobj == EV::e_system_cost);
#pragma endregion

#pragma region Electricity Network Constraints
	// C1, C2: number of generation units at each node	
	//int existP = 0;
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int i = 0; i < nPlt; i++)
	//	{
	//		// is plant type i part of initial plants types of node n:
	//		bool isInd = std::find(Enodes[n].Init_plt_ind.begin(), Enodes[n].Init_plt_ind.end(), i) != Enodes[n].Init_plt_ind.end();
	//		if (isInd)
	//		{
	//			// find the index of plant in the set of Init_plants of the node
	//			string tp = pltType2sym[i];
	//			int ind1 = std::find(Enodes[n].Init_plt_types.begin(), Enodes[n].Init_plt_types.end(), tp) - Enodes[n].Init_plt_types.begin();
	//			Model.addConstr(EV::Xop[n][i] == Enodes[n].Init_plt_count[ind1] - EV::Xdec[n][i] + EV::Xest[n][i]);
	//		}
	//		else
	//		{
	//			Model.addConstr(EV::Xop[n][i] == -EV::Xdec[n][i] + EV::Xest[n][i]);
	//		}
	//		//C2: maximum number of each plant type at each node (not necessary)
	//		//Model.addConstr(EV::Xop[n][i] <= Plants[i].Umax);
	//	}
	//}

	//C3, C4: production limit, ramping	
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Params::Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				//Model.addConstr(EV::prod[n][t][i] >= Plants[i].Pmin * EV::Xop[n][i]); since we don't consider unit commitment in this model
				Model.addConstr(EV::prod[n][t][i] <= Params::Plants[i].Pmax * EV::val_Xop[n][i]);

				//if (t > 0 && !Setting::heuristics1_active) 
				if (t > 0)
				{
					Model.addConstr(-Params::Plants[i].rampU * Params::Plants[i].Pmax * EV::val_Xop[n][i] <= EV::prod[n][t][i] - EV::prod[n][t - 1][i]);
					Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Params::Plants[i].rampU * Params::Plants[i].Pmax * EV::val_Xop[n][i]);
				}
			}
		}
	}

	// C 3.5: max yearly generation for existing plants
	for (int i = 0; i < nPlt; i++)
	{
		GRBLinExpr ex0(0);
		if (Params::Plants[i].is_exis != 1) { continue; }

		for (int n = 0; n < nEnode; n++)
		{
			for (int t = 0; t < Params::Te.size(); t++)
			{
				ex0 += Params::time_weight[t] * EV::prod[n][t][i];
			}
		}
		Model.addConstr(ex0 <= Params::Plants[i].max_yearly_gen);
	}




	// C5, C6: flow limit for electricity
	for (int br = 0; br < nBr; br++)
	{
		for (int t = 0; t < Params::Te.size(); t++)
		{
			if (Params::Branches[br].is_exist == 1)
			{
				Model.addConstr(EV::flowE[br][t] <= Params::Branches[br].maxFlow);
				Model.addConstr(-Params::Branches[br].maxFlow <= EV::flowE[br][t]);
			}
			else
			{
				Model.addConstr(EV::flowE[br][t] <= Params::Branches[br].maxFlow * EV::val_Ze[br]);
				Model.addConstr(-Params::Branches[br].maxFlow * EV::Ze[br] <= EV::flowE[br][t]);
			}
		}
	}

	// peak demand
	double max_dem = 0;
	for (int t = 0; t < Params::Te.size(); t++)
	{
		double demm = 0;
		for (int n = 0; n < nEnode; n++)
		{
			demm += Params::Enodes[n].demand[Params::Te[t]];
		}
		if (demm > max_dem)
		{
			max_dem = demm;
		}
	}
	//C7: power balance
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Params::Te.size(); t++)
		{
			GRBLinExpr exp_prod(0);
			for (int i = 0; i < nPlt; i++)
			{
				exp_prod += EV::prod[n][t][i];
			}
			GRBLinExpr exp_trans(0);
			for (int m : Params::Enodes[n].adj_buses)
			{
				// instead of defining a dictionray, one could traverse the entire branches...
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
				ex_store += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r];
			}
			double dem = Params::Enodes[n].demand[Params::Te[t]];
			//dem = dem/1.896;
			// ignore trans
			//Model.addConstr(exp_prod + ex_store + EV::curtE[n][t] == dem);
			//Model.addConstr(exp_prod + exp_trans +  EV::curtE[n][t] == dem);

			EV::PB[n][t] = Model.addConstr(exp_prod + exp_trans + ex_store + EV::curtE[n][t] == dem);
		}
	}
	// C8: flow equation
	int ebr = 0;
	for (int br = 0; br < Params::Branches.size(); br++)
	{
		int fb = Params::Branches[br].from_bus;
		int tb = Params::Branches[br].to_bus;
		for (int t = 0; t < Params::Te.size(); t++)
		{
			//if (Branches[br].is_exist == 1 && !Setting::heuristics1_active)
			if (Params::Branches[br].is_exist == 1)
			{
				Model.addConstr(EV::flowE[br][t] - Params::Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) <= 1e-2);
				Model.addConstr(-1e-2 <= EV::flowE[br][t] - Params::Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]));
			}
			//else if (!Setting::heuristics1_active)
			else
			{
				// Big M = Branches[br].suscep * 200
				Model.addConstr(EV::flowE[br][t] - Params::Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) <= 1e-2 + Params::Branches[br].suscep * 200 * (1 - EV::val_Ze[br]));
				Model.addConstr(-1e-2 - Params::Branches[br].suscep * 200 * (1 - EV::val_Ze[br]) <= EV::flowE[br][t] - Params::Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]));
			}
		}
	}

	// C9: phase angle (theta) limits. already applied in the definition of the variable
	for (int n = 0; n < nEnode; n++)
	{
		Model.addConstr(EV::theta[n][0] == 0);
		for (int t = 1; t < Params::Te.size(); t++)
		{
			Model.addConstr(EV::theta[n][t] <= pi);
			Model.addConstr(-pi <= EV::theta[n][t]);
		}
	}

	// C10: VRE production profile
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Params::Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Params::Plants[i].type == "solar" | Params::Plants[i].type == "solar-UPV")
				{
					Model.addConstr(EV::prod[n][t][i] <= Params::Plants[i].zonal_profile[Params::Te[t]][n] * Params::Plants[i].Pmax * EV::val_Xop[n][i]);
				}
				if (Params::Plants[i].type == "wind" || Params::Plants[i].type == "wind-new")
				{
					Model.addConstr(EV::prod[n][t][i] <= Params::Plants[i].zonal_profile[Params::Te[t]][n] * Params::Plants[i].Pmax * EV::val_Xop[n][i]);
				}
				if (Params::Plants[i].type == "wind_offshore" || Params::Plants[i].type == "wind-offshore-new")
				{
					Model.addConstr(EV::prod[n][t][i] <= Params::Plants[i].zonal_profile[Params::Te[t]][0] * Params::Plants[i].Pmax * EV::val_Xop[n][i]);
				}
			}
		}
	}


	// C11: demand curtainlment constraint
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Params::Te.size(); t++)
		{
			double dem = Params::Enodes[n].demand[Params::Te[t]];
			//dem = dem/1.896;
			Model.addConstr(EV::curtE[n][t] <= dem);
		}
	}

	//C12: RPS constraints
	GRBLinExpr exp3(0);
	double total_demand = 0;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Params::Te.size(); t++)
		{
			total_demand += Params::time_weight[t] * Params::Enodes[n].demand[Params::Te[t]];

			for (int i = 0; i < nPlt; i++)
			{
				//exp2 += time_weight[t] * Plants[i].emis_rate * EV::prod[n][t][i];
				if (Params::Plants[i].type == "solar" || Params::Plants[i].type == "wind" ||
					Params::Plants[i].type == "wind_offshore" || Params::Plants[i].type == "wind-offshore-new" ||
					Params::Plants[i].type == "hydro" || Params::Plants[i].type == "solar-UPV" ||
					Params::Plants[i].type == "wind-new" || Params::Plants[i].type == "hydro-new")
				{
					exp3 += Params::time_weight[t] * EV::prod[n][t][i]; // production from VRE
				}
			}
		}
	}
	Model.addConstr(exp3 >= Setting::RPS * total_demand);

	// C14,C15,C16 storage constraints
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		for (int r = 0; r < neSt; r++)
	//		{
	//			/*if (t == Te.size() - 1)
	//			{
	//				Model.addConstr(EV::eSlev[n][t][r] <= EV::eSlev[n][0][r]);
	//			}*/
	//			if (t == 0)
	//			{
	//				Model.addConstr(EV::eSlev[n][t][r] ==
	//					Estorage[r].eff_ch * EV::eSch[n][t][r] -
	//					(EV::eSdis[n][t][r] / Estorage[r].eff_disCh));
	//			}
	//			else
	//			{
	//				Model.addConstr(EV::eSlev[n][t][r] == EV::eSlev[n][t - 1][r] +
	//					Estorage[r].eff_ch * EV::eSch[n][t][r] -
	//					(EV::eSdis[n][t][r] / Estorage[r].eff_disCh));
	//			}
	//			Model.addConstr(EV::YeCD[n][r] >= EV::eSdis[n][t][r]);
	//			Model.addConstr(EV::YeCD[n][r] >= EV::eSch[n][t][r]);
	//			Model.addConstr(EV::YeLev[n][r] >= EV::eSlev[n][t][r]);
	//		}
	//	}
	//}

	//// start and ending of storage should be the same in case of using rep. days
	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int r = 0; r < neSt; r++)
	//	{
	//		for (int k = 0; k < Tg.size(); k++)
	//		{
	//			int str = k * 24;
	//			int end = (k + 1) * 24 - 1;
	//			Model.addConstr(EV::eSlev[n][str][r] == EV::eSlev[n][end][r]);
	//		}
	//	}
	//}


	// C17: if an storage established or not
	/*for (int n = 0; n < nEnode; n++)
	{
		for (int r = 0; r < neSt; r++)
		{
			Model.addConstr(EV::YeLev[n][r] <= 10e8 * EV::YeStr[n][r]);
			Model.addConstr(EV::YeCD[n][r] <= EV::YeLev[n][r]);
		}
	}*/



#pragma endregion

}


void Primal_subproblem()
{
	//auto start = chrono::high_resolution_clock::now();
#pragma region Fetch Data

	std::map<string, int> sym2pltType = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7},
		{"Ct",8},{"CC",9},{"CC-CCS",10},{"solar-UPV",11},{"wind-new",12},
		{"wind-offshore-new",13},{"hydro-new",14},{"nuclear-new",15} };
	std::map<int, string> pltType2sym = { {0,"ng"},{1,"dfo"},
{2,"solar"},{3,"wind"},{4,"wind_offshore"},{5,"hydro"},{6,"coal"},{7,"nuclear"} };
	/*vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<branch> Branches = Params::Branches;*/
	int nEnode = (int)Params::Enodes.size();
	int nPlt = (int)Params::Plants.size();
	int nBr = (int)Params::Branches.size();
	int neSt = (int)Params::Estorage.size();
	//vector<int> Tg = Params::Tg;
	//vector<int> Te = Params::Te;
	//vector<int> time_weight = Params::time_weight;
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
	//map<int, vector<int>> Le = Params::Le;
	//vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion

	auto start = chrono::high_resolution_clock::now();
	GRBEnv* envB = 0;
	envB = new GRBEnv();
	GRBModel ModelB = GRBModel(envB);
	GRBLinExpr exp_NGobj(0);
	GRBLinExpr exp_Eobj(0);
	Populate_EV(ModelB); Populate_GV(ModelB);
	Elec_Module_Primal_SP(ModelB, exp_Eobj);

	NG_Module(ModelB, exp_NGobj);

	// coupling constraints
	GRBLinExpr ex_xi(0);
	GRBLinExpr ex_NG_emis(0);
	GRBLinExpr ex_E_emis(0);
	//Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);

	double flowge = 0;
	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Params::Tg.size(); tau++)
		{
			for (int n : Params::Gnodes[k].adjE)
			{
				GRBLinExpr exp2(0);
				for (int t = tau * 24; t < (tau + 1) * 24; t++)
				{
					for (int i = 0; i < nPlt; i++)
					{
						if (Params::Plants[i].type == "ng" || Params::Plants[i].type == "CT" || Params::Plants[i].type == "CC" || Params::Plants[i].type == "CC-CCS")
						{

							exp2 += Params::time_weight[t] * Params::Plants[i].heat_rate * EV::prod[n][t][i];

						}
					}
				}
				ModelB.addConstr(Params::RepDaysCount[tau] * GV::val_flowGE[k][n][tau] - exp2 <= 1e-2);
				ModelB.addConstr(exp2 - Params::RepDaysCount[tau] * GV::val_flowGE[k][n][tau] <= 1e-2);
				ex_xi += Params::RepDaysCount[tau] * GV::val_flowGE[k][n][tau];
			}
		}
	}

	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Params::Tg.size(); tau++)
		{
			ex_NG_emis += Params::RepDaysCount[tau] * Params::NG_emis_rate * (GV::val_supply[k][tau]);
		}
	}
	ex_NG_emis -= Params::NG_emis_rate * CV::xi;

	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Params::Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Params::Plants[i].type == "ng" || Params::Plants[i].type == "CT" || Params::Plants[i].type == "CC" || Params::Plants[i].type == "CC-CCS")
				{
					ex_E_emis += Params::time_weight[t] * Params::Plants[i].emis_rate * Params::Plants[i].heat_rate * EV::prod[n][t][i];
				}
			}
		}
	}
	if (Setting::is_xi_given)
	{
		cout << "\n\n if xi given" << endl;
		ModelB.addConstr(ex_xi == Setting::xi_val * Setting::PGC);
		ModelB.addConstr(CV::xi == ex_xi);
	}
	else
	{
		cout << "else xi given" << endl;
		ModelB.addConstr(CV::xi == ex_xi);
	}

	/*if (Setting::Case == 1)
	{
		Model.addConstr(100 == CV::E_emis);
		Model.addConstr(100 == CV::NG_emis);
	}*/
	if (Setting::Case == 2)
	{
		ModelB.addConstr(ex_E_emis <= Setting::Emis_lim * Setting::PE);
		ModelB.addConstr(ex_E_emis == CV::E_emis);
		ModelB.addConstr(ex_NG_emis == CV::NG_emis);
	}
	if (Setting::Case == 3)
	{// the original model
		ModelB.addConstr(ex_E_emis + ex_NG_emis <= Setting::Emis_lim * Setting::PE);
		ModelB.addConstr(ex_E_emis == CV::E_emis);
		ModelB.addConstr(ex_NG_emis == CV::NG_emis);
	}

	ModelB.setObjective(exp_Eobj + exp_NGobj, GRB_MINIMIZE);

	ModelB.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	ModelB.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	//Model.set(GRB_IntParam_DualReductions, 0);
	ModelB.optimize();
	if (ModelB.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || ModelB.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize IM!!!" << endl;
		std::cout << ModelB.get(GRB_IntAttr_Status);
	}


	double obj_val = ModelB.get(GRB_DoubleAttr_ObjVal);
	double gap = 0;
	if (!Setting::relax_int_vars)
	{
		double gap = ModelB.get(GRB_DoubleAttr_MIPGap);
	}

	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "Elapsed time: " << Elapsed << endl;
	std::cout << "\t Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << ModelB.get(GRB_IntAttr_Status) << endl;
}




void Subproblem()
{
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
	GRBEnv* env = 0;
	env = new GRBEnv();
	GRBModel Model = GRBModel(env);
#pragma region Define Decision Variables
	SP::alpha = new GRBVar[nPlt];
	SP::beta = new GRBVar * *[nEnode];
	SP::gamma1 = new GRBVar * *[nEnode];
	SP::gamma2 = new GRBVar * *[nEnode];
	SP::delta11 = new GRBVar * [nBr];
	SP::delta12 = new GRBVar * [nBr];
	SP::delta21 = new GRBVar * [nBr];
	SP::delta22 = new GRBVar * [nBr];
	SP::theta = new GRBVar * [nEnode];
	SP::zeta11 = new GRBVar * [nBr];
	SP::zeta21 = new GRBVar * [nBr];
	SP::zeta12 = new GRBVar * [nBr];
	SP::zeta22 = new GRBVar * [nBr];
	SP::eta1 = new GRBVar * [nEnode];
	SP::eta2 = new GRBVar * [nEnode];
	SP::eta3 = new GRBVar[nEnode];
	SP::pi = new GRBVar * *[nEnode];
	SP::phi = new GRBVar * [nEnode];
	SP::rho = new GRBVar * *[nGnode];
	SP::omega = Model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	SP::tau = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);


	for (int i = 0; i < nPlt; i++)
	{
		SP::alpha[i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
	}

	for (int n = 0; n < nEnode; n++)
	{
		SP::beta[n] = new GRBVar * [Te.size()];
		SP::gamma1[n] = new GRBVar * [Te.size()]; SP::gamma2[n] = new GRBVar * [Te.size()];
		SP::theta[n] = new GRBVar[Te.size()];
		SP::eta1[n] = new GRBVar[Te.size()];
		SP::eta2[n] = new GRBVar[Te.size()];

		SP::eta3[n] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);

		SP::pi[n] = new GRBVar * [Te.size()];
		SP::phi[n] = new GRBVar[Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			SP::beta[n][t] = new GRBVar[nPlt];
			SP::gamma1[n][t] = new GRBVar[nPlt]; SP::gamma2[n][t] = new GRBVar[nPlt];

			// theta is a free dual variable
			SP::theta[n][t] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);

			SP::eta1[n][t] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::eta2[n][t] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
			SP::pi[n][t] = new GRBVar[nPlt];
			SP::phi[n][t] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
			for (int i = 0; i < nPlt; i++)
			{
				SP::beta[n][t][i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
				SP::gamma1[n][t][i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
				SP::gamma2[n][t][i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
				SP::pi[n][t][i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			}
		}
	}


	for (int b = 0; b < nBr; b++)
	{
		SP::delta11[b] = new GRBVar[Te.size()]; SP::delta12[b] = new GRBVar[Te.size()];
		SP::delta21[b] = new GRBVar[Te.size()]; SP::delta22[b] = new GRBVar[Te.size()];
		SP::zeta11[b] = new GRBVar[Te.size()]; SP::zeta12[b] = new GRBVar[Te.size()];
		SP::zeta22[b] = new GRBVar[Te.size()]; SP::zeta22[b] = new GRBVar[Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			SP::delta11[b][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			SP::delta12[b][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			SP::delta21[b][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			SP::delta22[b][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);

			SP::zeta11[b][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			SP::zeta12[b][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			SP::zeta21[b][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
			SP::zeta22[b][t] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
		}
	}

	for (int k = 0; k < nGnode; k++)
	{
		SP::rho[k] = new GRBVar * [nEnode];
		for (int n = 0; n < nEnode; n++)
		{
			SP::rho[k][n] = new GRBVar[Tg.size()];
			for (int tau = 0; tau < Tg.size(); tau++)
			{
				SP::rho[k][n][tau] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
			}
		}
	}
#pragma endregion



#pragma region Objective Function
	GRBLinExpr ex0(0);

	for (int i = 0; i < nPlt; i++) { ex0 += Plants[i].max_yearly_gen * SP::alpha[i]; }


	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				ex0 += Plants[i].Pmax * EV::val_Xop[n][i] * SP::beta[n][t][i];
				ex0 += Plants[i].rampU * EV::val_Xop[n][i] * SP::gamma1[n][t][i];
				ex0 += Plants[i].rampU * EV::val_Xop[n][i] * SP::gamma2[n][t][i];

				if (Plants[i].type == "solar" | Plants[i].type == "solar-UPV")
				{
					ex0 += Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::val_Xop[n][i] * SP::pi[n][t][i];
				}
				if (Plants[i].type == "wind" || Plants[i].type == "wind-new")
				{
					ex0 += Plants[i].zonal_profile[Te[t]][n] * Plants[i].Pmax * EV::val_Xop[n][i] * SP::pi[n][t][i];
				}
				if (Plants[i].type == "wind_offshore" || Plants[i].type == "wind-offshore-new")
				{
					ex0 += Plants[i].zonal_profile[Te[t]][0] * Plants[i].Pmax * EV::val_Xop[n][i] * SP::pi[n][t][i];
				}


			}
			double dem = Enodes[n].demand[Te[t]];
			double str = 0;;
			for (int r = 0; r < neSt; r++)
			{
				str += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r];
			}
			ex0 += (dem - str) * SP::theta[n][t];
			ex0 += pi * SP::eta1[n][t] + pi * SP::eta2[n][t];

			ex0 += dem * SP::phi[n][t];
			ex0 += Setting::RPS * dem * SP::omega;
		}
	}

	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (Branches[b].is_exist == 1)
			{
				ex0 += Branches[b].maxFlow * (SP::delta11[b][t] + SP::delta12[b][t]);

			}
			else
			{
				ex0 += Branches[b].maxFlow * EV::val_Ze[b] * (SP::delta21[b][t] + SP::delta22[b][t]);
				ex0 += Branches[b].suscep * 200 * (1 - EV::val_Ze[b]) * (SP::zeta11[b][t] + SP::zeta12[b][t]);
			}
		}
	}

	double rem = Setting::Emis_lim * Setting::PE;

	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			for (int n : Gnodes[k].adjE)
			{
				ex0 += GV::val_flowGE[k][n][tau] * SP::rho[k][n][tau];
				rem -= (RepDaysCount[tau] * Params::NG_emis_rate * (GV::val_supply[k][tau] - GV::val_flowGE[k][n][tau]));
			}
		}
	}
	ex0 += rem * SP::tau;

	Model.setObjective(ex0, GRB_MAXIMIZE);
#pragma endregion


#pragma region Constraints for prod variable
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "nuclear")
				{
					Model.addConstr(SP::alpha[i] + SP::beta[n][t][i] +
						SP::gamma1[n][t][i] - SP::gamma1[n][t + 1][i] -
						SP::gamma2[n][t][i] + SP::gamma2[n][t + 1][i] + SP::theta[n][t] +
						SP::pi[n][t][i] + SP::omega
						<= time_weight[t] * (Plants[i].var_cost + Params::nuclear_price));
				}
				if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
				{
					GRBLinExpr ex_adj(0);
					for (int k = 0; k < nGnode; k++)
					{
						for (int n2 : Gnodes[k].adjE)
						{
							if (n2 == n)
							{
								int tau = t / 24;
								ex_adj += time_weight[t] * Plants[i].heat_rate * SP::rho[k][n2][tau];
							}
						}
					}

					Model.addConstr(SP::alpha[i] + SP::beta[n][t][i] +
						SP::gamma1[n][t][i] - SP::gamma1[n][t + 1][i] -
						SP::gamma2[n][t][i] + SP::gamma2[n][t + 1][i] + SP::theta[n][t] +
						SP::pi[n][t][i] + SP::omega + ex_adj +
						time_weight[t] * Plants[i].emis_rate * Plants[i].heat_rate * SP::tau
						<= time_weight[t] * (Plants[i].var_cost));
				}
				else
				{
					{
						Model.addConstr(SP::alpha[i] + SP::beta[n][t][i] +
							SP::gamma1[n][t][i] - SP::gamma1[n][t + 1][i] -
							SP::gamma2[n][t][i] + SP::gamma2[n][t + 1][i] + SP::theta[n][t] +
							SP::pi[n][t][i] + SP::omega
							<= time_weight[t] * (Plants[i].var_cost));
					}
				}
			}
		}
	}
#pragma endregion

#pragma region Constraints for flowE variable (flowE is a free variable)
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			GRBLinExpr ex_fe(0);
			if (Branches[b].is_exist == 1)
			{
				ex_fe += SP::delta11 - SP::delta12;
				ex_fe += SP::zeta11 - SP::zeta12;
			}
			else
			{
				ex_fe += SP::delta21 - SP::delta22;
				ex_fe += SP::zeta21 - SP::zeta22;
			}
			if (Branches[b].from_bus > Branches[b].to_bus)
			{
				ex_fe -= SP::theta[Branches[b].from_bus][t];
			}
			else
			{
				ex_fe += SP::theta[Branches[b].from_bus][t];
			}

			Model.addConstr(ex_fe == 0);

		}
	}
#pragma endregion


#pragma region Constratints for phase angle, load shedding, 
	// load shedding
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			Model.addConstr(SP::theta[n][t] + SP::phi[n][t] <= time_weight[t] * E_curt_cost);
		}
	}

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
							exp_line += -Branches[l2].suscep * (SP::zeta11[l2][t] - SP::zeta12[l2][t]);
						}
						else
						{
							exp_line += -Branches[l2].suscep * (-SP::zeta11[l2][t] + SP::zeta12[l2][t]);

						}
					}
					else
					{
						if (n > m)
						{
							exp_line += -Branches[l2].suscep * (SP::zeta21[l2][t] - SP::zeta22[l2][t]);
						}
						else
						{
							exp_line += -Branches[l2].suscep * (-SP::zeta21[l2][t] + SP::zeta22[l2][t]);
						}
					}
				}
			}

			if (t == 0)
			{
				Model.addConstr(exp_line + SP::eta1[n][t] - SP::eta2[n][t] + SP::eta3[n] <= 0);
			}
			else
			{
				Model.addConstr(exp_line + SP::eta1[n][t] - SP::eta2[n][t] <= 0);
			}
		}
	}
#pragma endregion

	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	//Model.set(GRB_IntParam_DualReductions, 0);
	Model.optimize();
	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize IM!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
	}


	double obj_val = Model.get(GRB_DoubleAttr_ObjVal);
	double gap = Model.get(GRB_DoubleAttr_MIPGap);

	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "Elapsed time: " << Elapsed << endl;
	std::cout << "\t Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;

}




