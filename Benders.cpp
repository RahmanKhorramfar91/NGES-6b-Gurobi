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
	double pi = 3.141592;
	int nGnode = (int)Gnodes.size();
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
	GRBEnv* env = 0;
	env = new GRBEnv();
	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_NGobj(0);
	GRBLinExpr exp_Eobj(0);
	Populate_EV(Model); Populate_GV(Model);
	Elec_Module_Primal_SP(Model, exp_Eobj);

	NG_Module(Model, exp_NGobj);

	// coupling constraints
	GRBLinExpr ex_xi(0);
	GRBLinExpr ex_NG_emis(0);
	GRBLinExpr ex_E_emis(0);
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);

	double flowge = 0;
	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			for (int n : Gnodes[k].adjE)
			{
				GRBLinExpr exp2(0);
				for (int t = tau * 24; t < (tau + 1) * 24; t++)
				{
					for (int i = 0; i < nPlt; i++)
					{
						if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
						{
							
								exp2 += time_weight[t] * Plants[i].heat_rate * EV::prod[n][t][i];
							
						}
					}
				}
				Model.addConstr(RepDaysCount[tau] * GV::val_flowGE[k][n][tau] - exp2 <= 1e-2);
				Model.addConstr(exp2 - RepDaysCount[tau] * GV::val_flowGE[k][n][tau] <= 1e-2);
				ex_xi += RepDaysCount[tau] * GV::val_flowGE[k][n][tau];
			}
		}
	}

	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			ex_NG_emis += RepDaysCount[tau] * Params::NG_emis_rate * (GV::val_supply[k][tau]);
		}
	}
	ex_NG_emis -= Params::NG_emis_rate * CV::xi;

	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
				{
					ex_E_emis += time_weight[t] * Plants[i].emis_rate * Plants[i].heat_rate * EV::prod[n][t][i];
				}
			}
		}
	}
	if (Setting::is_xi_given)
	{
		cout << "\n\n if xi given" << endl;
		Model.addConstr(ex_xi == Setting::xi_val * Setting::PGC);
		Model.addConstr(CV::xi == ex_xi);
	}
	else
	{
		cout << "else xi given" << endl;
		Model.addConstr(CV::xi == ex_xi);
	}

	/*if (Setting::Case == 1)
	{
		Model.addConstr(100 == CV::E_emis);
		Model.addConstr(100 == CV::NG_emis);
	}*/
	if (Setting::Case == 2)
	{
		Model.addConstr(ex_E_emis <= Setting::Emis_lim * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}
	if (Setting::Case == 3)
	{// the original model
		Model.addConstr(ex_E_emis + ex_NG_emis <= Setting::Emis_lim * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}

	Model.setObjective(exp_Eobj + exp_NGobj, GRB_MINIMIZE);

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
	double gap = 0;
	if (!Setting::relax_int_vars)
	{
		double gap = Model.get(GRB_DoubleAttr_MIPGap);
	}

	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "Elapsed time: " << Elapsed << endl;
	std::cout << "\t Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;
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




