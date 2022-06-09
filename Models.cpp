#include"Models_Funcs.h"

double feas_sol(double& elec_LB, double& elec_UB, double& ng_obj, double& feas_gap)
{
	auto start = chrono::high_resolution_clock::now();

#pragma region Fetch Data
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	int nEnode = (int)Enodes.size();
	int nPlt = (int)Plants.size();
	vector<int> Te = Params::Te;
#pragma endregion

	Setting::heuristics1_active = true; Setting::relax_int_vars = true;
	Setting::warm_start_active = false;

#pragma region Feasible solution for electricity network
	// solve heuristic to get a LB

	GRBEnv* env = 0;
	env = new GRBEnv();
	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_Eobj(0);
	Populate_EV_SP(Model); Populate_GV(Model);
	Elec_Module(Model, exp_Eobj);
	Model.setObjective(exp_Eobj, GRB_MINIMIZE);
	GRBLinExpr ex_xi(0);
	GRBLinExpr ex_NG_emis(0);
	GRBLinExpr ex_E_emis(0);
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);
	// Coupling 1
	Model.addConstr(CV::xi == ex_xi);
	Model.addConstr(ex_xi == Setting::xi_val * Setting::PGC);

	// Coupling 2
	Model.addConstr(ex_E_emis <= Setting::Emis_lim * Setting::PE);
	Model.addConstr(ex_E_emis == CV::E_emis);

	Model.set(GRB_DoubleParam_TimeLimit, 3600);
	Model.set(GRB_DoubleParam_MIPGap, 0.01);
	Model.optimize();
	double gapp = Model.get(GRB_DoubleAttr_MIPGap);
	elec_LB = Model.get(GRB_DoubleAttr_ObjBound);
	/*double gg0 = Model.get(GRB_DoubleAttr_MinBound);
	double gg2 = Model.get(GRB_DoubleParam_BestBdStop);
	double gg3 = Model.get(GRB_DoubleAttr_BoundVio);
	double gg = Model.get(GRB_DoubleAttr_MinBound);*/
	Get_EV_vals(Model);
	Setting::heuristics1_active = false; Setting::relax_int_vars = false;
	env->~GRBEnv();
	Model.~GRBModel();


	// set Xop variables and solve the model to get an UB
	GRBEnv* env2 = 0;
	env2 = new GRBEnv();
	GRBModel Model2 = GRBModel(env2);
	exp_Eobj.clear();
	Populate_EV_SP(Model2); Populate_GV(Model2);
	Elec_Module(Model2, exp_Eobj);
	Model2.setObjective(exp_Eobj, GRB_MINIMIZE);
	ex_xi.clear();
	ex_NG_emis.clear();
	ex_E_emis.clear();
	Coupling_Constraints(Model2, ex_xi, ex_NG_emis, ex_E_emis);
	// Coupling 1
	Model2.addConstr(CV::xi == ex_xi);
	Model2.addConstr(ex_xi == Setting::xi_val * Setting::PGC);

	// Coupling 2
	Model2.addConstr(ex_E_emis <= Setting::Emis_lim * Setting::PE);
	Model2.addConstr(ex_E_emis == CV::E_emis);
	//Model.addConstr(CV::NG_emis);

	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			Model2.addConstr(EV::Xop[n][i] == std::round(EV::val_Xop[n][i]));
			//Model2.addConstr(EV::Xest[n][i] == EV::val_Xest[n][i]);
			//Model2.addConstr(EV::Xdec[n][i] == EV::val_Xdec[n][i]);			
		}
	}
	for (int b = 0; b < (int)Params::Branches.size(); b++)
	{
		Model2.addConstr(EV::Ze[b] == std::round(EV::val_Ze[b]));
	}
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < (int)Params::Estorage.size(); i++)
		{
			Model2.addConstr(EV::YeStr[n][i] == std::round(EV::val_YeStr[n][i]));
		}
	}

	Model2.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model2.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model2.optimize();

	/*if (Model2.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
		Setting::warm_start_active = false;
		return 0;
	}*/


	elec_UB = Model2.get(GRB_DoubleAttr_ObjVal);
	feas_gap = Model2.get(GRB_DoubleAttr_MIPGap);
	CV::used_emis_cap = CV::E_emis.get(GRB_DoubleAttr_X);
	Get_EV_vals(Model2);
	Model2.~GRBModel();
	env2->~GRBEnv();
#pragma endregion

#pragma region Feasible solutin for NG network
	GRBEnv* env3 = 0;
	env3 = new GRBEnv();
	GRBModel Model3 = GRBModel(env3);
	GRBLinExpr exp_NGobj3(0);
	Populate_GV(Model3);
	Populate_EV_SP(Model3);
	NG_Module(Model3, exp_NGobj3);

	// coupling constraints
	Setting::DGSP_active = true;
	GRBLinExpr ex_xi3(0); GRBLinExpr ex_NG_emis3(0); GRBLinExpr ex_E_emis3(0);
	Coupling_Constraints(Model3, ex_xi3, ex_NG_emis3, ex_E_emis3);
	Setting::DGSP_active = false;//back to default
	// Coupling 1: xi must be given
	Model3.addConstr(CV::xi == ex_xi3);
	//Model3.addConstr(ex_xi3 == Setting::xi_val*Setting::PGC);


	// Coupling 2
	double rhs = Setting::Emis_lim * Setting::PE - CV::used_emis_cap;
	rhs = std::max(rhs, 0 + 0.001);
	Model3.addConstr(ex_NG_emis3 <= rhs);
	Model3.addConstr(ex_NG_emis3 == CV::NG_emis);

	Model3.setObjective(exp_NGobj3, GRB_MINIMIZE);

	Model3.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model3.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model3.optimize();

	ng_obj = Model3.get(GRB_DoubleAttr_ObjVal);
	double rt = elec_UB + ng_obj;
	Get_GV_vals(Model3);
	//Model3.~GRBModel();
#pragma endregion

	Setting::warm_start_active = true;

	return 0;
}

void Electricy_Network_Model(GRBEnv* env)
{
	auto start = chrono::high_resolution_clock::now();

#pragma region apply heuristic to generate a warm-start solution
	double elec_UB = 0;
	if (Setting::warm_start_active)
	{
		double elec_LB = 0; double elec_UB = 0; double ng_obj; double feas_gap = 0;
		elec_UB = feas_sol(elec_LB, elec_UB, ng_obj, feas_gap);
	}
#pragma endregion

	/*GRBEnv* env = 0;
	env = new GRBEnv();*/
	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_Eobj(0);
	Populate_EV_SP(Model);	Populate_GV(Model);
	Elec_Module(Model, exp_Eobj);
	Model.setObjective(exp_Eobj, GRB_MINIMIZE);
	CV::xi = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	CV::NG_emis = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	CV::E_emis = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	Model.addConstr(CV::xi == 1);
	//Model.addConstr(CV::E_emis == 1);
	//Model.addConstr(CV::NG_emis == 1);
	// coupling constraints
	GRBLinExpr ex_xi(0);
	GRBLinExpr ex_NG_emis(0);
	GRBLinExpr ex_E_emis(0);
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);


	//if (Setting::Case == 1)
	{
		Model.addConstr(ex_E_emis <= Setting::Emis_lim * Setting::PE);
		Model.addConstr(ex_E_emis == CV::E_emis);
		Model.addConstr(ex_NG_emis == CV::NG_emis);
	}
	

#pragma region Solve the model

	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model.optimize();

	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize Electricity network model!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
		throw(-1);
	}

	double obj_val = Model.get(GRB_DoubleAttr_ObjVal);
	double gap = 0;
	if (!Setting::relax_int_vars)
	{
		gap = Model.get(GRB_DoubleAttr_MIPGap);

	}
	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "Elapsed time: " << Elapsed << endl;
	std::cout << "\t Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;
#pragma endregion



	Get_EV_vals(Model);

}

void NG_Network_Model(GRBEnv* env)
{
	auto start = chrono::high_resolution_clock::now();

	/*GRBEnv* env = 0;
	env = new GRBEnv();*/
	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_NGobj(0);

	Populate_GV(Model);
	NG_Module(Model, exp_NGobj);
	Model.setObjective(exp_NGobj, GRB_MINIMIZE);

	CV::xi = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	CV::NG_emis = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	CV::E_emis = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	Model.addConstr(CV::xi == 100);
	Model.addConstr(CV::E_emis == 100);
	Model.addConstr(CV::NG_emis == 100);

#pragma region Solve the model
	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model.optimize();

	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
		throw(-1);
	}

	double obj_val = Model.get(GRB_DoubleAttr_ObjVal);
	double gap = 0;
	if (!Setting::relax_int_vars)
	{
		gap = Model.get(GRB_DoubleAttr_MIPGap);
	}
	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "Elapsed time: " << Elapsed << endl;
	std::cout << "\t Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;
#pragma endregion
	Get_GV_vals(Model);
}

double Integrated_Model(GRBEnv* env)
{
	auto start = chrono::high_resolution_clock::now();
	 
#pragma region apply heuristic to generate a warm-start solution
	double elec_UB = 0;
	if (Setting::warm_start_active)
	{
		double elec_LB = 0; double elec_UB = 0; double ng_obj; double feas_gap = 0;
		elec_UB = feas_sol(elec_LB, elec_UB, ng_obj, feas_gap);
	}
#pragma endregion
	/*GRBEnv* env = 0;
	env = new GRBEnv();*/
	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_NGobj(0);
	GRBLinExpr exp_Eobj(0);
	Populate_EV_SP(Model);
	Elec_Module(Model, exp_Eobj);

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


	if (Setting::xi_LB_obj)
	{
		Model.setObjective(ex_xi, GRB_MINIMIZE);
	}
	else if (Setting::xi_UB_obj)
	{
		Model.setObjective(ex_xi, GRB_MAXIMIZE);
	}
	else
	{
		Model.setObjective(exp_Eobj + exp_NGobj, GRB_MINIMIZE);
		if (Setting::warm_start_active)
		{
			Model.addConstr(exp_Eobj + exp_NGobj <= elec_UB);
		}
	}

#pragma region Solve the model
	//Model.set(GRB_IntParam_OutputFlag, 0);
	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	//Model.set(GRB_IntParam_DualReductions, 0);

	Model.optimize();
	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize IM!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
		return -1;
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
	std::cout << "\t xi value = " << CV::xi.get(GRB_DoubleAttr_X) << endl;
	std::cout << "\t NG emission = " << CV::NG_emis.get(GRB_DoubleAttr_X) << endl;
	std::cout << "\t E emission = " << CV::E_emis.get(GRB_DoubleAttr_X) << endl;
#pragma endregion

	EV::MIP_gap = 0;
	GV::MIP_gap = 0;
	if (!Setting::relax_int_vars)
	{
		EV::MIP_gap = Model.get(GRB_DoubleAttr_MIPGap);
		GV::MIP_gap = Model.get(GRB_DoubleAttr_MIPGap);
	}
	Get_EV_vals(Model);
	Get_GV_vals(Model);
	return obj_val;
}

double DESP(GRBEnv* env)
{
	auto start = chrono::high_resolution_clock::now();

	//GRBEnv* env = 0;
	//env = new GRBEnv();
	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_Eobj(0);
	Populate_EV_SP(Model); Populate_GV(Model);
	Elec_Module(Model, exp_Eobj);

	// coupling constraints
	GRBLinExpr ex_xi(0); GRBLinExpr ex_NG_emis(0); GRBLinExpr ex_E_emis(0);
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);

	// Coupling 1
	Model.addConstr(CV::xi == ex_xi);
	Model.addConstr(ex_xi == Setting::xi_val * Setting::PGC);

	// Coupling 2
	Model.addConstr(ex_E_emis <= Setting::Emis_lim * Setting::PE);
	//Model.addConstr(ex_NG_emis <= 2.6e6);
	Model.addConstr(ex_E_emis == CV::E_emis);
	//Model.addConstr(CV::NG_emis);


	Model.setObjective(exp_Eobj, GRB_MINIMIZE);


#pragma region Solve the model

	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model.optimize();
	double gap = Model.get(GRB_DoubleAttr_MIPGap);
	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize DESP!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
		return -1;
	}

	double obj_val = Model.get(GRB_DoubleAttr_ObjVal);
	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "Elapsed time: " << Elapsed << endl;
	std::cout << "\t DESP Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;
#pragma endregion
	CV::used_emis_cap = CV::E_emis.get(GRB_DoubleAttr_X);

	Get_EV_vals(Model);

	return obj_val;
}

double DGSP(GRBEnv* env)
{
	auto start = chrono::high_resolution_clock::now();

	/*GRBEnv* env = 0;
	env = new GRBEnv();*/
	GRBModel Model = GRBModel(env);
	GRBLinExpr exp_NGobj(0);
	Populate_GV(Model);
	Populate_EV_SP(Model);
	NG_Module(Model, exp_NGobj);

	// coupling constraints
	Setting::DGSP_active = true;
	GRBLinExpr ex_xi(0); GRBLinExpr ex_NG_emis(0); GRBLinExpr ex_E_emis(0);
	Coupling_Constraints(Model, ex_xi, ex_NG_emis, ex_E_emis);
	Setting::DGSP_active = false;//back to default
	// Coupling 1: xi must be given
	Model.addConstr(CV::xi == ex_xi);
	//Model.addConstr(ex_xi == Setting::xi_val*Setting::PGC);


	// Coupling 2
	double rhs = Setting::Emis_lim * Setting::PE - CV::used_emis_cap;
	rhs = std::max(rhs, 0 + 0.001);
	Model.addConstr(ex_NG_emis <= rhs);
	Model.addConstr(ex_NG_emis == CV::NG_emis);

	Model.setObjective(exp_NGobj, GRB_MINIMIZE);

#pragma region Solve the model
	Model.set(GRB_DoubleParam_TimeLimit, Setting::CPU_limit);
	Model.set(GRB_DoubleParam_MIPGap, Setting::cplex_gap);
	Model.optimize();

	if (Model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || Model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
	{
		std::cout << "Failed to optimize DGSP!!!" << endl;
		std::cout << Model.get(GRB_IntAttr_Status);
		return -1;
	}

	double obj_val = Model.get(GRB_DoubleAttr_ObjVal);
	double gap = Model.get(GRB_DoubleAttr_MIPGap);
	auto end = chrono::high_resolution_clock::now();
	double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
	std::cout << "Elapsed time: " << Elapsed << endl;
	std::cout << "\t DGSP Obj Value:" << obj_val << endl;
	std::cout << "\t Gap: " << gap << " Status:" << Model.get(GRB_IntAttr_Status) << endl;
	std::cout << "\t xi value = " << CV::xi.get(GRB_DoubleAttr_X) << endl;
	std::cout << "\t NG emission = " << CV::NG_emis.get(GRB_DoubleAttr_X) << endl;
	std::cout << "\t E emission = " << CV::E_emis.get(GRB_DoubleAttr_X) << endl;
#pragma endregion
	Get_GV_vals(Model);
	return obj_val;
}

