#include"Models_Funcs.h";

// validate Benders SP by solving the primal Subproblem
void SP_flow_Upper()
{
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
			exp_Eobj += flowE_abs1[b][t] ;
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
	for (int i = 0; i < nPlt; i++)
	{
		GRBLinExpr ex0(0);
		if (Plants[i].is_exis != 1) { continue; }

		for (int n = 0; n < nEnode; n++)
		{
			for (int t = 0; t < Te.size(); t++)
			{
				ex0 += time_weight[t] * EV::prod[n][t][i];
			}
		}
		Model.addConstr(ex0 <= Plants[i].max_yearly_gen);
	}
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
		Model.addConstr(EV::theta[n][0] == 0);
		for (int t = 1; t < Te.size(); t++)
		{
			Model.addConstr(-EV::theta[n][t] >= -pi);
			Model.addConstr(EV::theta[n][t] >= -pi);
		}
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
void Primal_subproblem(vector<SP>& Cuts)
{
	double eps = 0.02;
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

#pragma region Objective Function
	//Model.addConstr(exp_Eobj == 11);
	exp_Eobj = 11;
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
				SP::d_delta11[br][t] = Model.addConstr(-EV::flowE[br][t] >= -Branches[br].maxFlow);
				SP::d_delta12[br][t] = Model.addConstr(EV::flowE[br][t] >= -Branches[br].maxFlow);
			}
			else
			{
				double rhs = (Branches[br].maxFlow * EV::val_Ze[br]);
				SP::d_delta21[br][t] = Model.addConstr(-EV::flowE[br][t] >= -rhs);
				SP::d_delta22[br][t] = Model.addConstr(EV::flowE[br][t] >= -rhs);
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

			SP::d_theta1[n][t] = Model.addConstr(exp_trans >= rhs);
			SP::d_theta2[n][t] = Model.addConstr(-exp_trans >= -rhs);
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
		Model.addConstr(EV::theta[n][0] == 0);
		for (int t = 1; t < Te.size(); t++)
		{
			SP::d_eta1[n][t] = Model.addConstr(-EV::theta[n][t] >= -pi);
			SP::d_eta2[n][t] = Model.addConstr(EV::theta[n][t] >= -pi);
		}
	}

	// C10: VRE production profile


	//// C11: demand curtainlment constraint

	////C12: RPS constraints
	//// C14,C15,C16 storage constraints


#pragma endregion


	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);

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
	std::cout << "\n\n Primal SP Elapsed time: " << Elapsed << endl;
	std::cout << "\t Primal SP Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;


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
	SP::dual_val_delta11 = new double* [nBr];
	SP::dual_val_delta12 = new double* [nBr];
	SP::dual_val_delta21 = new double* [nBr];
	SP::dual_val_delta22 = new double* [nBr];
	SP::dual_val_theta1 = new double* [nEnode];
	SP::dual_val_theta2 = new double* [nEnode];
	SP::dual_val_zeta11 = new double* [nBr];
	SP::dual_val_zeta12 = new double* [nBr];
	SP::dual_val_zeta21 = new double* [nBr];
	SP::dual_val_zeta22 = new double* [nBr];
	SP::dual_val_eta1 = new double* [nEnode];
	SP::dual_val_eta2 = new double* [nEnode];

	for (int n = 0; n < nEnode; n++)
	{
		SP::dual_val_theta1[n] = new double[Te.size()]();
		SP::dual_val_theta2[n] = new double[Te.size()]();

		SP::dual_val_eta1[n] = new double[Te.size()]();
		SP::dual_val_eta2[n] = new double[Te.size()]();
		double ave_price = 0;
		for (int t = 0; t < Te.size(); t++)
		{
			SP::dual_val_theta1[n][t] = std::max(0.0, SP::d_theta1[n][t].get(GRB_DoubleAttr_Pi));
			SP::dual_val_theta2[n][t] = std::max(0.0, SP::d_theta2[n][t].get(GRB_DoubleAttr_Pi));
			if (t > 0)
			{
				SP::dual_val_eta1[n][t] = std::max(0.0, SP::d_eta1[n][t].get(GRB_DoubleAttr_Pi));
				SP::dual_val_eta2[n][t] = std::max(0.0, SP::d_eta2[n][t].get(GRB_DoubleAttr_Pi));
			}
		}
	}

	for (int b = 0; b < nBr; b++)
	{
		SP::dual_val_delta11[b] = new double[Te.size()];
		SP::dual_val_delta12[b] = new double[Te.size()];
		SP::dual_val_delta21[b] = new double[Te.size()];
		SP::dual_val_delta22[b] = new double[Te.size()];

		SP::dual_val_zeta11[b] = new double[Te.size()];
		SP::dual_val_zeta12[b] = new double[Te.size()];
		SP::dual_val_zeta21[b] = new double[Te.size()];
		SP::dual_val_zeta22[b] = new double[Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			if (Branches[b].is_exist == 1)
			{
				SP::dual_val_delta11[b][t] = std::max(0.0, SP::d_delta11[b][t].get(GRB_DoubleAttr_Pi));
				SP::dual_val_delta12[b][t] = std::max(0.0, SP::d_delta12[b][t].get(GRB_DoubleAttr_Pi));
				SP::dual_val_zeta21[b][t] = std::max(0.0, SP::d_zeta21[b][t].get(GRB_DoubleAttr_Pi));
				SP::dual_val_zeta22[b][t] = std::max(0.0, SP::d_zeta22[b][t].get(GRB_DoubleAttr_Pi));
			}
			else
			{
				SP::dual_val_delta21[b][t] = std::max(0.0, SP::d_delta21[b][t].get(GRB_DoubleAttr_Pi));
				SP::dual_val_delta22[b][t] = std::max(0.0, SP::d_delta22[b][t].get(GRB_DoubleAttr_Pi));


				SP::dual_val_zeta11[b][t] = std::max(0.0, SP::d_zeta11[b][t].get(GRB_DoubleAttr_Pi));
				SP::dual_val_zeta12[b][t] = std::max(0.0, SP::d_zeta12[b][t].get(GRB_DoubleAttr_Pi));
			}
		}
	}

	std::cout << "dual_val_theta1: " << endl;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (SP::dual_val_theta1[n][t] > 0.0)
			{
				std::cout << SP::dual_val_theta1[n][t] << endl;
			}
		}
	}
	std::cout << "dual_val_theta2: " << endl;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (SP::dual_val_theta2[n][t] > 0.0)
			{
				std::cout << SP::dual_val_theta2[n][t] << endl;
			}
		}
	}


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

	//double rem = Setting::Emis_lim * Setting::PE;

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

	//rem = Setting::Emis_lim * Setting::PE;

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
	SP sp_new;
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
	}

	std::cout << "dual_val_eta1: " << endl;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (sp_new.dual_val_eta1[n][t] > 0.0)
			{
				std::cout << sp_new.dual_val_eta1[n][t] << endl;
			}
		}
	}
	std::cout << "dual_val_eta2: " << endl;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (sp_new.dual_val_eta2[n][t] > 0.0)
			{
				std::cout << sp_new.dual_val_eta2[n][t] << endl;
			}
		}
	}
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_theta[n][t] = EV::theta[n][t].get(GRB_DoubleAttr_X);
			if (EV::val_theta[n][t] > 0.0)
			{
				std::cout << EV::theta[n][t].get(GRB_DoubleAttr_X) << endl;
			}
		}
	}

	std::cout << "val_flowE" << endl;
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_flowE[b][t] = EV::flowE[b][t].get(GRB_DoubleAttr_X);
			if (EV::val_flowE[b][t] > 0)
			{
				//std::cout << "flowE[" << b << "][" << t << "] = " << EV::val_flowE[b][t] << endl;
			}
		}
	}


	for (int br = 0; br < 10; br++)
	{
		int fb = Branches[br].from_bus;
		int tb = Branches[br].to_bus;
		for (int t = 0; t < 3; t++)
		{
			if (Branches[br].is_exist == 1)
			{
				//Model.addConstr(EV::flowE[br][t] >= Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) - eps);
				//Model.addConstr(-EV::flowE[br][t] >= -Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) - eps);
				//std::cout << "theta[from]: " << EV::val_theta[fb][t] << "\t theta[to]: " << EV::val_theta[tb][t] << endl;
				//std::cout << "flowE: " << EV::val_flowE[br][t] << "\t rhs: " << Branches[br].suscep * (EV::val_theta[tb][t] - EV::val_theta[fb][t]) << endl;

			}
		}
	}
#pragma endregion



	Cuts.push_back(sp_new);
	Setting::Approach_1_active = true; Setting::Approach_2_active = false;
}



void Benders_Decomposition()
{
	vector<SP> Cuts;
	//Primal_subproblem(Cuts);
	//Dual_Subproblem();
	Master_Problem(Cuts);

	Primal_subproblem(Cuts);
}

void Master_Problem(vector<SP> Cuts)
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
	GRBEnv* env = 0;
	env = new GRBEnv();
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
		std::cout << "\n\n if xi given" << endl;
		Model.addConstr(ex_xi == Setting::xi_val * Setting::PGC);
		Model.addConstr(CV::xi == ex_xi);
	}
	else
	{
		std::cout << "else xi given" << endl;
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




	// define the contribution from subproblem
	GRBVar* Chi;
	GRBLinExpr ex_cut(0);
	Chi = Model.addVars(2, GRB_CONTINUOUS);
	Model.addConstr(Chi[0] >= ex_cut);
	//Chi[0] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);

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

			if (Plants[i].type == "nuclear-new")
			{
				Model.addConstr(EV::Xest[n][i] == 0);
			}

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

				// emission cost (no longer needed)
				/*if (Plants[i].type != "ng" && Plants[i].type != "CT" && Plants[i].type != "CC" && Plants[i].type != "CC-CCS")
				{
					ex_emis += time_weight[t] * Plants[i].emis_cost * Plants[i].emis_rate * EV::prod[n][t][i];
				}*/
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
	exp_Eobj = ex_est + ex_decom + ex_fix + ex_emis + ex_var + ex_thermal_fuel + ex_shedd + ex_trans + ex_elec_str;

	Model.addConstr(EV::est_cost == ex_est);
	Model.addConstr(EV::decom_cost == ex_decom);
	Model.addConstr(EV::fixed_cost == ex_fix);
	Model.addConstr(EV::var_cost == ex_var);
	Model.addConstr(EV::thermal_fuel_cost == ex_thermal_fuel);
	Model.addConstr(EV::dfo_coal_emis_cost == ex_emis);
	Model.addConstr(EV::shedding_cost == ex_shedd);
	Model.addConstr(EV::elec_storage_cost == ex_elec_str);

	Model.addConstr(exp_Eobj == EV::e_system_cost);
	Model.setObjective(exp_Eobj + exp_NGobj, GRB_MINIMIZE);
#pragma endregion


#pragma region MP Electricity local constraints
	// C1, C2: number of generation units at each node	
	int existP = 0;
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
			//C2: maximum number of each plant type at each node (not necessary)
			//Model.addConstr(EV::Xop[n][i] <= Plants[i].Umax);
		}
	}

	//C3, C4: production limit, ramping	
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				//Model.addConstr(EV::prod[n][t][i] >= Plants[i].Pmin * EV::Xop[n][i]); since we don't consider unit commitment in this model
				Model.addConstr(EV::prod[n][t][i] <= Plants[i].Pmax * EV::Xop[n][i]);

				//if (t > 0 && !Setting::heuristics1_active) 
				if (t > 0)
				{
					//Model.addConstr(Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i] >= -EV::prod[n][t][i] + EV::prod[n][t - 1][i]);
					Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
					Model.addConstr(-EV::prod[n][t][i] + EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);

				}
			}
		}
	}

	// C 3.5: max yearly generation for existing plants
	for (int i = 0; i < nPlt; i++)
	{
		GRBLinExpr ex0(0);
		if (Plants[i].is_exis != 1) { continue; }

		for (int n = 0; n < nEnode; n++)
		{
			for (int t = 0; t < Te.size(); t++)
			{
				ex0 += time_weight[t] * EV::prod[n][t][i];
			}
		}
		Model.addConstr(ex0 <= Plants[i].max_yearly_gen);
	}

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


#pragma region Add cuts
	//GRBLinExpr ex_cut(0);
	ex_cut = 0;
	for (int c = 0; c < Cuts.size(); c++)
	{
		//for (int i = 0; i < nPlt; i++)
		//{
		//	if (Plants[i].is_exis != 1) { continue; } // only applies to existing plants
		//	ex_cut -= Plants[i].max_yearly_gen * Cuts[c].dual_val_alpha[i];
		//}
		double lhs = 0;
		for (int n = 0; n < nEnode; n++)
		{
			for (int t = 0; t < Te.size(); t++)
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

				ex_cut += (dem - str - exp_prod - EV::curtE[n][t]) * Cuts[c].dual_val_theta1[n][t];
				ex_cut -= (dem - str - exp_prod - EV::curtE[n][t]) * Cuts[c].dual_val_theta2[n][t];

				ex_cut -= pi * Cuts[c].dual_val_eta1[n][t];
				ex_cut -= pi * Cuts[c].dual_val_eta2[n][t];
				lhs += pi * Cuts[c].dual_val_eta1[n][t];
				lhs += pi * Cuts[c].dual_val_eta2[n][t];
				if (Cuts[c].dual_val_eta2[n][t] > 0)
				{
					std::cout << Cuts[c].dual_val_eta2[n][t] << endl;
				}
			}
		}

		for (int b = 0; b < nBr; b++)
		{
			for (int t = 0; t < Te.size(); t++)
			{
				if (Branches[b].is_exist == 1)
				{
					ex_cut -= Branches[b].maxFlow * Cuts[c].dual_val_delta11[b][t];
					ex_cut -= Branches[b].maxFlow * Cuts[c].dual_val_delta12[b][t];
					lhs += Branches[b].maxFlow * Cuts[c].dual_val_delta11[b][t];
					lhs += Branches[b].maxFlow * Cuts[c].dual_val_delta12[b][t];
					if (Cuts[c].dual_val_delta12[b][t] > 0)
					{
						std::cout << Cuts[c].dual_val_delta12[b][t] << endl;
					}
				}
				else
				{
					double Big_M = std::abs(Branches[b].suscep * 24);
					ex_cut -= Branches[b].maxFlow * EV::Ze[b] * Cuts[c].dual_val_delta21[b][t];
					ex_cut -= Branches[b].maxFlow * EV::Ze[b] * Cuts[c].dual_val_delta22[b][t];
					ex_cut += Big_M * (1 - EV::Ze[b]) * Cuts[c].dual_val_zeta11[b][t];
					ex_cut -= Big_M * (1 - EV::Ze[b]) * Cuts[c].dual_val_zeta12[b][t];
				}
			}
		}

		//double rem = Setting::Emis_lim * Setting::PE;
		//GRBLinExpr ex_em(0);
		/*ex_em = Setting::Emis_lim * Setting::PE;
		for (int k = 0; k < nGnode; k++)
		{
			for (int tau = 0; tau < Tg.size(); tau++)
			{
				GRBLinExpr fge(0);
				for (int n : Gnodes[k].adjE)
				{
					ex_cut += GV::flowGE[k][n][tau] * Cuts[c].dual_val_rho1[k][n][tau];
					ex_cut -= GV::flowGE[k][n][tau] * Cuts[c].dual_val_rho2[k][n][tau];
					fge += GV::flowGE[k][n][tau];
				}
				ex_em -= RepDaysCount[tau] * Params::NG_emis_rate * (GV::supply[k][tau] - fge);
			}
		}

		ex_cut -= ex_em * Cuts[c].dual_val_tau;*/

	}
	Model.addConstr(Chi[0] >= ex_cut);
	//Model.addConstr(Chi >= SP::SP_Primal_obj);
#pragma endregion



#pragma region Solve MP
	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model.set(GRB_IntParam_OutputFlag, 0);
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
	std::cout << "\t MP Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;
	std::cout << "\t Chi value = " << Chi[0].get(GRB_DoubleAttr_X) << endl;
	std::cout << "\t NG emission = " << CV::NG_emis.get(GRB_DoubleAttr_X) << endl;
	std::cout << "\t E emission = " << CV::E_emis.get(GRB_DoubleAttr_X) << endl;
#pragma endregion

	Get_EV_vals(Model);
	Get_GV_vals(Model);
	Print_Results(0, obj_val);

#pragma region Get the values of MP variables
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
			EV::val_num_est[i] += EV::Xest[n][i].get(GRB_DoubleAttr_X);
			EV::val_num_decom[i] += EV::Xdec[n][i].get(GRB_DoubleAttr_X);
			EV::val_Xop[n][i] = std::round(EV::Xop[n][i].get(GRB_DoubleAttr_X));
			EV::val_Xest[n][i] = std::round(EV::Xest[n][i].get(GRB_DoubleAttr_X));
			EV::val_Xdec[n][i] = std::round(EV::Xdec[n][i].get(GRB_DoubleAttr_X));
			if (EV::val_Xop[n][i] > 0)
			{
				std::cout << n << "\t" << i << "\t" << EV::val_Xop[n][i] << endl;
			}
		}
	}
	EV::val_Ze = new double[nBr]();
	for (int b = 0; b < nBr; b++)
	{
		EV::val_num_est_trans += EV::Ze[b].get(GRB_DoubleAttr_X);
		EV::val_Ze[b] = std::round(EV::Ze[b].get(GRB_DoubleAttr_X));
	}

	EV::val_eSch = new double** [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_eSch[n] = new double* [Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_eSch[n][t] = new double[neSt]();
			for (int i = 0; i < neSt; i++)
			{
				EV::val_eSch[n][t][i] = EV::eSch[n][t][i].get(GRB_DoubleAttr_X);
				//if (eSchS[n][t][i] > 10e-3)
				//{
				//	//	std::cout << "prod[" << n << "][" << t << "][" << i << "] = " << prodS[n][t][i] << endl;
				//	fid << "eS_ch[" << n << "][" << t << "][" << i << "] = " << eSchS[n][t][i] << endl;
				//}
			}
		}
	}
	EV::val_eSdis = new double** [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_eSdis[n] = new double* [Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_eSdis[n][t] = new double[neSt]();
			for (int i = 0; i < neSt; i++)
			{
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
	}
#pragma endregion


#pragma region Check Chi val
	double Chi_val = 0;

	for (int c = 0; c < Cuts.size(); c++)
	{
		//for (int i = 0; i < nPlt; i++)
		//{
		//	if (Plants[i].is_exis != 1) { continue; } // only applies to existing plants
		//	Chi_val -= Plants[i].max_yearly_gen * Cuts[c].dual_val_alpha[i];
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
				Chi_val += (dem - str - exp_prod - EV::val_curtE[n][t]) * Cuts[c].dual_val_theta1[n][t];
				Chi_val -= (dem - str - exp_prod - EV::val_curtE[n][t]) * Cuts[c].dual_val_theta2[n][t];
				Chi_val -= pi * SP::dual_val_eta1[n][t];
				Chi_val -= pi * SP::dual_val_eta2[n][t];
			}
		}

		for (int b = 0; b < nBr; b++)
		{
			for (int t = 0; t < Te.size(); t++)
			{
				if (Branches[b].is_exist == 1)
				{
					Chi_val -= Branches[b].maxFlow * Cuts[c].dual_val_delta11[b][t];
					Chi_val -= Branches[b].maxFlow * Cuts[c].dual_val_delta12[b][t];
				}
				else
				{
					double Big_M = std::abs(Branches[b].suscep * 24);
					Chi_val -= Branches[b].maxFlow * EV::val_Ze[b] * Cuts[c].dual_val_delta21[b][t];
					Chi_val -= Branches[b].maxFlow * EV::val_Ze[b] * Cuts[c].dual_val_delta22[b][t];
					Chi_val += Big_M * (1 - EV::val_Ze[b]) * Cuts[c].dual_val_zeta11[b][t];
					Chi_val -= Big_M * (1 - EV::val_Ze[b]) * Cuts[c].dual_val_zeta12[b][t];
				}
			}
		}
	}
	std::cout << "Chi rhs: " << Chi_val << endl;
#pragma endregion


}







void Dual_Subproblem()
{
	double eps = 0.4;
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
	SP::eta3 = new GRBVar[nEnode];
	/*SP::pi = new GRBVar * *[nEnode];
	SP::phi = new GRBVar * [nEnode];
	SP::rho = new GRBVar * *[nGnode];
	SP::omega = Model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	SP::tau = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);*/


	/*for (int i = 0; i < nPlt; i++)
	{
		SP::alpha[i] = Model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS);
	}*/

	for (int n = 0; n < nEnode; n++)
	{
		/*SP::beta[n] = new GRBVar * [Te.size()];
		SP::gamma1[n] = new GRBVar * [Te.size()]; SP::gamma2[n] = new GRBVar * [Te.size()];*/
		SP::theta1[n] = new GRBVar[Te.size()];
		SP::theta2[n] = new GRBVar[Te.size()];
		SP::eta1[n] = new GRBVar[Te.size()];
		SP::eta2[n] = new GRBVar[Te.size()];

		// free variable
		SP::eta3[n] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);

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
			/*Model.addConstr(SP::delta21[b][t] == 0);
			Model.addConstr(SP::delta22[b][t] == 0);
			Model.addConstr(SP::delta11[b][t] == 0);
			Model.addConstr(SP::delta12[b][t] == 0);

			Model.addConstr(SP::zeta11[b][t] == 0);
			Model.addConstr(SP::zeta12[b][t] == 0);
			Model.addConstr(SP::zeta21[b][t] == 0);
			Model.addConstr(SP::zeta22[b][t] == 0);*/
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
			ex0 -= (dem - str - exp_prod - EV::val_curtE[n][t]) * SP::theta1[n][t];
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
				ex0 += Big_M * (1 - EV::val_Ze[b]) * SP::zeta11[b][t];
				ex0 -= Big_M * (1 - EV::val_Ze[b]) * SP::zeta12[b][t];
			}
		}
	}



	//ex0 = 11;
	//ex0 = ex0 + SP::tau;

	Model.setObjective(ex0, GRB_MAXIMIZE);
	//Model.addConstr(ex0 == 9.68247e+08);
	Model.addConstr(ex0 <= 0);
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
	//for (int b = 0; b < nBr; b++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		double ex_fe = 0;
	//		if (Branches[b].is_exist == 1)
	//		{
	//			ex_fe += SP::dual_val_delta11[b][t] - SP::dual_val_delta12[b][t];
	//			ex_fe += SP::dual_val_zeta21[b][t] - SP::dual_val_zeta22[b][t];
	//		}
	//		else
	//		{
	//			ex_fe += SP::dual_val_delta21[b][t] - SP::dual_val_delta22[b][t];
	//			ex_fe += SP::dual_val_zeta11[b][t] - SP::dual_val_zeta12[b][t];
	//		}
	//		/*if (Branches[b].from_bus > Branches[b].to_bus)
	//		{
	//			ex_fe -= SP::theta[Branches[b].from_bus][t];
	//		}
	//		else
	//		{
	//			ex_fe += SP::theta[Branches[b].from_bus][t];
	//		}*/
	//		ex_fe += (SP::dual_val_theta[Branches[b].from_bus][t] - SP::dual_val_theta[Branches[b].to_bus][t]);
	//		//std::cout << "from bus: " << Branches[b].from_bus << " \t to bus: " << Branches[b].to_bus << endl;
	//		
	//		if (ex_fe<-0.001 || ex_fe>0.001)
	//		{
	//			int gg = 0;
	//		}
	//	}
	//}


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
				ex_fe += SP::zeta11[b][t];
				ex_fe -= SP::zeta12[b][t];
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
							exp_line -= Branches[l2].suscep * (SP::zeta11[l2][t] - SP::zeta12[l2][t]);
						}
						else
						{
							exp_line += Branches[l2].suscep * (SP::zeta11[l2][t] - SP::zeta12[l2][t]);

						}
					}
					else
					{
						if (n > m)
						{
							exp_line -= Branches[l2].suscep * (SP::zeta21[l2][t] - SP::zeta22[l2][t]);
						}
						else
						{
							exp_line += Branches[l2].suscep * (-SP::zeta21[l2][t] + SP::zeta22[l2][t]);
						}
					}
				}
			}

			if (t == 0)
			{
				Model.addConstr(exp_line - SP::eta1[n][t] + SP::eta2[n][t] + SP::eta3[n] == 0);
			}
			else
			{
				Model.addConstr(exp_line + SP::eta1[n][t] - SP::eta2[n][t] == 0);
			}
		}
	}
#pragma endregion

	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	//Model.set(GRB_IntParam_DualReductions, 0);
	Model.set(GRB_IntParam_OutputFlag, 0);
	Model.optimize();
	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize IM!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
	}
	/*if (~Model.get(GRB_IntAttr_Status)==GRB_OPTIMAL)
	{
		std::cout << "Failed to optimize IM!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
	}*/

	double obj_val = Model.get(GRB_DoubleAttr_ObjVal);

	//double gap = Model.get(GRB_DoubleAttr_MIPGap);

	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "\n\n\n\n Dual Elapsed time: " << Elapsed << endl;
	std::cout << "\t Dual Obj Value:" << obj_val << endl;
	std::cout << " Dual Status:" << Model.get(GRB_IntAttr_Status) << endl;


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
	for (int n = 0; n < nEnode; n++)
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
	}

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
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (t > 10)
			{
				//break;
			}
			double dz = SP::delta11[b][t].get(GRB_DoubleAttr_X);
			if (dz < -0.001 || dz>0.001)
			{
				std::cout << "delta11[" << b << "][" << t << "] = " << dz << endl;
			}
			dz = SP::delta12[b][t].get(GRB_DoubleAttr_X);
			if (dz < -0.001 || dz>0.001)
			{
				std::cout << "delta12[" << b << "][" << t << "] = " << dz << endl;
			}

			dz = SP::delta21[b][t].get(GRB_DoubleAttr_X);
			if (dz < -0.001 || dz>0.001)
			{
				std::cout << "delta21[" << b << "][" << t << "] = " << dz << endl;
			}

			dz = SP::delta22[b][t].get(GRB_DoubleAttr_X);
			if (dz < -0.001 || dz>0.001)
			{
				std::cout << "delta22[" << b << "][" << t << "] = " << dz << endl;
			}
		}
	}

	////zeta
	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (t > 10)
			{
				//break;
			}
			double dz = SP::zeta11[b][t].get(GRB_DoubleAttr_X);
			if (dz < -0.001 || dz>0.001)
			{
				std::cout << "zeta11[" << b << "][" << t << "] = " << dz << endl;
			}
			dz = SP::zeta12[b][t].get(GRB_DoubleAttr_X);
			if (dz < -0.001 || dz>0.001)
			{
				std::cout << "zeta12[" << b << "][" << t << "] = " << dz << endl;
			}
			dz = SP::zeta21[b][t].get(GRB_DoubleAttr_X);
			if (dz < -0.001 || dz>0.001)
			{
				std::cout << "zeta21[" << b << "][" << t << "] = " << dz << endl;
			}
			dz = SP::zeta22[b][t].get(GRB_DoubleAttr_X);
			if (dz < -0.001 || dz>0.001)
			{
				std::cout << "zeta22[" << b << "][" << t << "] = " << dz << endl;
			}
		}
	}


	//eta
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (t > 10)
			{
				//break;
			}
			double eta = SP::eta1[n][t].get(GRB_DoubleAttr_X);
			if (eta < -0.001 || eta>0.001)
			{
				std::cout << "eta1[" << n << "][" << t << "] = " << eta << endl;
			}

			eta = SP::eta2[n][t].get(GRB_DoubleAttr_X);
			if (eta < -0.001 || eta>0.001)
			{
				std::cout << "eta1[" << n << "][" << t << "] = " << eta << endl;
			}

			eta = SP::eta3[n].get(GRB_DoubleAttr_X);
			if (eta < -0.001 || eta>0.001)
			{
				std::cout << "eta1[" << n << "][" << t << "] = " << eta << endl;
			}
		}
	}

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


	double obj_theta = 0; double obj_other = 0;
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{

			double dem = Enodes[n].demand[Te[t]];
			double str = 0;
			for (int r = 0; r < neSt; r++)
			{
				str += EV::val_eSdis[n][t][r] - EV::val_eSch[n][t][r];
			}
			//std::cout << "dem-str: " << dem - str << endl;
			obj_theta += (dem - str) * SP::theta1[n][t].get(GRB_DoubleAttr_X); // commented
			obj_other += pi * SP::eta1[n][t].get(GRB_DoubleAttr_X) + pi * SP::eta2[n][t].get(GRB_DoubleAttr_X);
		}
	}

	for (int b = 0; b < nBr; b++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (b == 27 && t == 16)
			{
				std::cout << EV::val_Ze[b] << endl;
			}
			if (Branches[b].is_exist == 1)
			{
				obj_other += Branches[b].maxFlow * (SP::delta11[b][t].get(GRB_DoubleAttr_X) + SP::delta12[b][t].get(GRB_DoubleAttr_X));
			}
			else
			{
				obj_other += Branches[b].maxFlow * EV::val_Ze[b] * (SP::delta21[b][t].get(GRB_DoubleAttr_X) + SP::delta22[b][t].get(GRB_DoubleAttr_X));
				obj_other += std::abs(Branches[b].suscep) * 24 * (1 - EV::val_Ze[b]) * (SP::zeta11[b][t].get(GRB_DoubleAttr_X) + SP::zeta12[b][t].get(GRB_DoubleAttr_X));
			}
		}
	}




	std::cout << "dual obj_theta value: " << obj_theta << endl;
	std::cout << "dual obj_other value: " << obj_other << endl;
	std::cout << "total dual obj value: " << obj_theta + obj_other << endl;
#pragma endregion

}




