#include"Models_Funcs.h"

#pragma region EV struct
GRBVar** EV::Xest; // integer (continues) for plants
GRBVar** EV::Xdec; // integer (continues) for plants
GRBVar** EV::YeCD; // continuous: charge/discharge capacity
GRBVar** EV::YeLev; // continuous: charge/discharge level
GRBVar** EV::YeStr;
GRBVar* EV::Ze;
GRBVar** EV::theta; // continuous phase angle
GRBVar** EV::curtE; // continuous curtailment variable
GRBVar*** EV::prod;// continuous
GRBVar*** EV::eSch;// power charge to storage 
GRBVar*** EV::eSdis;// power discharge to storage
GRBVar*** EV::eSlev;// power level at storage
GRBVar** EV::Xop; // integer (continues)
GRBVar** EV::flowE; // unlike the paper, flowE subscripts are "ntm" here
GRBVar EV::est_cost;
GRBVar EV::decom_cost;
GRBVar EV::fixed_cost;
GRBVar EV::var_cost;
GRBVar EV::thermal_fuel_cost;
GRBVar EV::shedding_cost;
GRBVar EV::elec_storage_cost;
GRBVar EV::Emit_var;
GRBVar EV::e_system_cost;
GRBVar EV::dfo_coal_emis_cost;

GRBConstr** EV::PB;

double*** EV::val_prod;
double*** EV::val_sCh;
double*** EV::val_sDis;
double EV::val_est_cost;
double EV::val_decom_cost;
double EV::val_fixed_cost;
double EV::val_var_cost;
double EV::val_thermal_fuel_cost;
double EV::val_shedding_cost;
double EV::val_elec_storage_cost;
double EV::val_Emit_var;
double EV::val_e_system_cost;
double EV::val_num_storage; // number of storage facilities established
double EV::val_storage_lev; // value of toal storage capacity (MWh)
double EV::val_storage_cap; // value of toal storage capacity (MW)
double* EV::val_num_est = new double[Params::Plants.size()]();// number of established plants
double* EV::val_num_decom = new double[Params::Plants.size()](); // number of plants decommissioned
double  EV::val_total_curt; // total load shedding
double  EV::val_num_est_trans; // number of new transmission lines established
double  EV::val_total_flow; // total flow in the E network
double* EV::val_total_prod = new double[Params::Plants.size()](); // total production	
double EV::val_dfo_coal_emis_cost;
double** EV::val_Xop;
double** EV::val_Xest;
double** EV::val_Xdec;
double* EV::val_Ze;
double EV::MIP_gap;
double** EV::val_YeStr;
double** EV::val_curtE;
double*** EV::val_eSlev;
double** EV::val_PB;
#pragma endregion

#pragma region NG struct
GRBVar* GV::Xstr;
GRBVar* GV::Xvpr;
GRBVar** GV::Sstr;
GRBVar** GV::Svpr;
GRBVar** GV::Sliq;
GRBVar** GV::supply;
GRBVar** GV::curtG;
GRBVar** GV::curtRNG;
GRBVar** GV::flowGG;
GRBVar*** GV::flowGE;
GRBVar*** GV::flowGL;
GRBVar*** GV::flowVG;
GRBVar* GV::Zg;
GRBVar GV::strInv_cost;
GRBVar GV::pipe_cost;
GRBVar GV::gShedd_cost;
GRBVar GV::rngShedd_cost;
GRBVar GV::gStrFOM_cost;
GRBVar GV::NG_import_cost;
GRBVar GV::NG_system_cost;
double*** GV::val_flowGE;
double GV::val_strInv_cost;
double GV::val_pipe_cost;
double GV::val_ngShedd_cost;
double GV::val_rngShedd_cost;
double GV::val_gStrFOM_cost;
double GV::val_NG_import_cost;
double GV::val_NG_system_cost;
double* GV::val_supply;
double GV::val_ng_curt;
double GV::val_rng_curt;
int GV::val_num_est_pipe;
double GV::val_total_flowGG;
double GV::val_total_flowGE;
double GV::val_total_flowGL;
double GV::val_total_flowVG;
double* GV::val_storage;
double* GV::val_vapor;
double* GV::val_Xstr;
double* GV::val_Zg;
double GV::MIP_gap;
#pragma endregion


#pragma region Coupling struct (CV)
GRBVar CV::xi;
GRBVar CV::NG_emis;
GRBVar CV::E_emis;
double CV::used_emis_cap;
double CV::val_xi; // f
double CV::val_NG_emis;
double CV::val_E_emis;
#pragma endregion

void  Populate_EV(GRBModel& Model)
{

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

#pragma region Electricity network DVs
	EV::Xest = new GRBVar * [nEnode]; // integer (continues) for plants
	EV::Xdec = new GRBVar * [nEnode]; // integer (continues) for plants
	EV::YeCD = new GRBVar * [nEnode]; // continuous: charge/discharge capacity
	EV::YeLev = new GRBVar * [nEnode]; // continuous: charge/discharge level
	EV::YeStr = new GRBVar * [nEnode]; // (binary) if a storage is established
	EV::PB = new GRBConstr * [nEnode];
	if (Setting::relax_int_vars)
	{
		EV::Ze = Model.addVars(nBr, GRB_CONTINUOUS);
	}
	else
	{
		EV::Ze = Model.addVars(nBr, GRB_BINARY);
	}

	EV::theta = new GRBVar * [nEnode]; // continuous phase angle
	EV::curtE = new GRBVar * [nEnode]; // continuous curtailment variable

	EV::prod = new GRBVar * *[nEnode]; // continuous
	EV::eSch = new GRBVar * *[nEnode];// power charge to storage 
	EV::eSdis = new GRBVar * *[nEnode];// power discharge to storage
	EV::eSlev = new GRBVar * *[nEnode];// power level at storage

	EV::Xop = new GRBVar * [nEnode];// integer (continues)
	EV::flowE = new GRBVar * [nBr]; // unlike the paper, flowE subscripts are "ntm" here
	for (int b = 0; b < nBr; b++)
	{
		EV::flowE[b] = Model.addVars(Te.size());
		for (int t = 0; t < Te.size(); t++)
		{
			EV::flowE[b][t] = Model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			//EV::flowE[b][t] = Model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		//Model.addVars()
	}
	for (int n = 0; n < nEnode; n++)
	{
		EV::PB[n] = new GRBConstr[Te.size()];
		if (Setting::relax_int_vars)
		{
			EV::Xest[n] = Model.addVars(nPlt, GRB_CONTINUOUS);
			EV::Xdec[n] = Model.addVars(nPlt, GRB_CONTINUOUS);
			EV::Xop[n] = Model.addVars(nPlt, GRB_CONTINUOUS);
		}
		else
		{
			EV::Xest[n] = Model.addVars(nPlt);
			EV::Xop[n] = Model.addVars(nPlt);
			EV::Xdec[n] = Model.addVars(nPlt);
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "solar-UPV" || Plants[i].type == "wind-new")
				{
					EV::Xest[n][i] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
					EV::Xdec[n][i] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
					EV::Xop[n][i] = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
				}
				else
				{
					EV::Xest[n][i] = Model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
					EV::Xdec[n][i] = Model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
					EV::Xop[n][i] = Model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
				}
			}
		}


		//Ze[n] = GRBVarArray(env, Enodes[n].adj_buses.size(), 0, IloInfinity, ILOBOOL);
		EV::theta[n] = Model.addVars((int)Te.size(), GRB_CONTINUOUS);
		EV::curtE[n] = Model.addVars((int)Te.size(), GRB_CONTINUOUS);
		EV::YeCD[n] = Model.addVars(neSt, GRB_CONTINUOUS);
		EV::YeLev[n] = Model.addVars(neSt, GRB_CONTINUOUS);
		if (Setting::relax_int_vars)
		{
			EV::YeStr[n] = Model.addVars(neSt, GRB_CONTINUOUS);

		}
		else
		{
			EV::YeStr[n] = Model.addVars(neSt, GRB_BINARY);
		}

		EV::prod[n] = new GRBVar * [Te.size()];
		EV::eSch[n] = new GRBVar * [Te.size()];
		EV::eSdis[n] = new GRBVar * [Te.size()];
		EV::eSlev[n] = new GRBVar * [Te.size()];

		for (int t = 0; t < Te.size(); t++)
		{
			EV::prod[n][t] = Model.addVars(nPlt, GRB_CONTINUOUS);
			EV::eSch[n][t] = Model.addVars(neSt, GRB_CONTINUOUS);
			EV::eSdis[n][t] = Model.addVars(neSt, GRB_CONTINUOUS);
			EV::eSlev[n][t] = Model.addVars(neSt, GRB_CONTINUOUS);
		}
	}
	EV::est_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::decom_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::fixed_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::var_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::thermal_fuel_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::dfo_coal_emis_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::shedding_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::elec_storage_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::Emit_var = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	EV::e_system_cost = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
#pragma endregion
}

void Populate_GV(GRBModel& Model)
{

#pragma region Fetch Data

	vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<exist_gSVL> Exist_SVL = Params::Exist_SVL;
	int nSVL = (int)Exist_SVL.size();
	vector<SVL> SVLs = Params::SVLs;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Enodes.size();
	int nPlt = (int)Plants.size();
	int nBr = (int)Branches.size();
	int neSt = (int)Estorage.size();
	int nPipe = (int)PipeLines.size();
	vector<int> Tg = Params::Tg;
	//vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	double pi = 3.1415926535897;
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
	map<int, vector<int>> Lg = Params::Lg;
	vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion

#pragma region NG network DVs

	// NG vars

	if (Setting::relax_int_vars)
	{
		GV::Xstr = Model.addVars(nSVL, GRB_CONTINUOUS);
		GV::Zg = Model.addVars(nPipe, GRB_CONTINUOUS);
	}
	else
	{
		GV::Xstr = Model.addVars(nSVL, GRB_BINARY);
		GV::Zg = Model.addVars(nPipe, GRB_BINARY);
	}

	GV::Xvpr = Model.addVars(nSVL, GRB_CONTINUOUS);

	GV::Sstr = new GRBVar * [nSVL];
	GV::Svpr = new GRBVar * [nSVL];
	GV::Sliq = new GRBVar * [nSVL];
	GV::supply = new GRBVar * [nGnode];
	GV::curtG = new GRBVar * [nGnode];
	GV::curtRNG = new GRBVar * [nGnode];
	GV::flowGG = new GRBVar * [nPipe];
	GV::flowGE = new GRBVar * *[nGnode];
	GV::flowGL = new GRBVar * *[nGnode];
	GV::flowVG = new GRBVar * *[nSVL];


	GV::strInv_cost = Model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	GV::pipe_cost = Model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	GV::gShedd_cost = Model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	GV::rngShedd_cost = Model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	GV::gStrFOM_cost = Model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	GV::NG_import_cost = Model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	GV::NG_system_cost = Model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);


	for (int i = 0; i < nPipe; i++)
	{
		// uni-directional
		GV::flowGG[i] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
	}
	for (int j = 0; j < nSVL; j++)
	{
		GV::Sstr[j] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
		GV::Svpr[j] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
		GV::Sliq[j] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
		GV::flowVG[j] = new GRBVar * [nGnode];
		for (int k = 0; k < nGnode; k++)
		{
			GV::flowVG[j][k] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
		}
	}

	for (int k = 0; k < nGnode; k++)
	{
		GV::supply[k] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
		GV::curtG[k] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
		GV::curtRNG[k] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
		GV::flowGE[k] = new GRBVar * [nEnode];
		GV::flowGL[k] = new GRBVar * [nSVL];
		for (int j = 0; j < nSVL; j++)
		{
			GV::flowGL[k][j] = Model.addVars(Tg.size(), GRB_CONTINUOUS);

		}
		for (int n = 0; n < nEnode; n++)
		{
			GV::flowGE[k][n] = Model.addVars(Tg.size(), GRB_CONTINUOUS);
		}
	}
#pragma endregion

}


void Elec_Module(GRBModel& Model, GRBLinExpr& exp_Eobj)
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


#pragma region Set some variables
	// 1) existing types can not be established because there are new equivalent types
	// 2) new types cannot be decommissioned
	// 3) no new wind-offshore on some location (not yet implemented)
	// 4) no new nuclear
	// 5) turn off df0/coal
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			//// do not allow establishment off-shore wind and hydro plants
			//Model.addConstr(EV::Xdec[n][5] == 0);
			if (Plants[i].type == "dfo" || Plants[i].type == "coal")
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
	
	/*GRBLinExpr exp0(0);
	for (int n = 0; n < nEnode; n++)
	{
		for (int i = 0; i < nPlt; i++)
		{
			if (Plants[i].type == "dfo" || Plants[i].type == "coal")
			{
				Model.addConstr(EV::Xop[n][i] == 0);
				continue;
			}
			if (!Plants[i].is_exis == 1)
			{
				Model.addConstr(EV::Xop[n][i] == 0);
			}
			else
			{
				Model.addConstr(EV::Xdec[n][i] == 0);
				Model.addConstr(EV::Xest[n][i] == 0);
			}
		}
		exp0 += EV::prod[n][0][0];
	}
	Model.addConstr(exp0 >= 0);*/
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
			if (Plants[i].type == "dfo" || Plants[i].type == "coal") { continue; }
			double s1 = (double)std::pow(1.0 / (1 + WACC), Plants[i].lifetime);
			double capco = WACC / (1 - s1) * Plants[i].Reg_coeffs_per_state[n];
			ex_est += capco * Plants[i].capex * EV::Xest[n][i];
			ex_decom += Plants[i].decom_cost * EV::Xdec[n][i];
		}

		// fixed cost (annual, so no iteration over time)
		for (int i = 0; i < nPlt; i++)
		{
			if (Plants[i].type == "dfo" || Plants[i].type == "coal") { continue; }
			ex_fix += Plants[i].fix_cost * EV::Xop[n][i];
		}
		// var+fuel costs of plants
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "dfo" || Plants[i].type == "coal") { continue; }
				// var cost
				ex_var += time_weight[t] * Plants[i].var_cost * EV::prod[n][t][i];

				// fuel price to be updated later (dollar per thousand cubic feet=MMBTu)
				// NG fuel is handle in the NG network for case 2 and 3
				if (Plants[i].type == "dfo")
				{
					ex_thermal_fuel += time_weight[t] * dfo_pric * Plants[i].heat_rate * EV::prod[n][t][i];

				}
				else if (Plants[i].type == "coal")
				{
					ex_thermal_fuel += time_weight[t] * coal_price * Plants[i].heat_rate * EV::prod[n][t][i];
				}
				else if (Plants[i].type == "nuclear" || Plants[i].type == "nuclear-new")
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
#pragma endregion

#pragma region Electricity Network Constraints
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
					Model.addConstr(-Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i] <= EV::prod[n][t][i] - EV::prod[n][t - 1][i]);
					Model.addConstr(EV::prod[n][t][i] - EV::prod[n][t - 1][i] <= Plants[i].rampU * Plants[i].Pmax * EV::Xop[n][i]);
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
	for (int br = 0; br < nBr; br++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			if (Branches[br].is_exist == 1)
			{
				Model.addConstr(EV::flowE[br][t] <= Branches[br].maxFlow);
				Model.addConstr(-Branches[br].maxFlow <= EV::flowE[br][t]);
			}
			else
			{
				Model.addConstr(EV::flowE[br][t] <= Branches[br].maxFlow * EV::Ze[br]);
				Model.addConstr(-Branches[br].maxFlow * EV::Ze[br] <= EV::flowE[br][t]);
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
				int key1 = n * 200 + m;
				// check if this key exist, if not, the other order exisit
				if (Le.count(key1) == 0)
				{
					key1 = m * 200 + n;
				}
				//each Le can contains multiple lines
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
			dem = dem;
			// ignore trans
			//Model.addConstr(exp_prod + ex_store + EV::curtE[n][t] == dem);
			//Model.addConstr(exp_prod + exp_trans +  EV::curtE[n][t] == dem);

			EV::PB[n][t] = Model.addConstr(exp_prod + exp_trans + ex_store + EV::curtE[n][t] == dem);
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
				Model.addConstr(EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) <= 1e-2);
				Model.addConstr(-1e-2 <= EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]));
			}
			//else if (!Setting::heuristics1_active)
			else
			{
				// Big M = Branches[br].suscep * 200
				Model.addConstr(EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]) <= 1e-2 + Branches[br].suscep * 200 * (1 - EV::Ze[br]));
				Model.addConstr(-1e-2 - Branches[br].suscep * 200 * (1 - EV::Ze[br]) <= EV::flowE[br][t] - Branches[br].suscep * (EV::theta[tb][t] - EV::theta[fb][t]));
			}
		}
	}

	// C9: phase angle (theta) limits. already applied in the definition of the variable
	for (int n = 0; n < nEnode; n++)
	{
		Model.addConstr(EV::theta[n][0] == 0);
		for (int t = 1; t < Te.size(); t++)
		{
			Model.addConstr(EV::theta[n][t] <= pi);
			Model.addConstr(-pi <= EV::theta[n][t]);
		}
	}

	// C10: VRE production profile
	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "solar" | Plants[i].type == "solar-UPV")
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
			dem = dem;
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

#pragma region add the warm start solution
	if (Setting::warm_start_active && !Setting::heuristics1_active)
	{
		for (int n = 0; n < nEnode; n++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				//Model.set(GRB_DoubleAttr_Start, EV::Xop[n],EV::val_Xopt[n],nPlt);
				EV::Xop[n][i].set(GRB_DoubleAttr_Start, EV::val_Xop[n][i]);
				EV::Xest[n][i].set(GRB_DoubleAttr_Start, EV::val_Xest[n][i]);
				EV::Xdec[n][i].set(GRB_DoubleAttr_Start, EV::val_Xdec[n][i]);
				for (int t = 0; t < Te.size(); t++)
				{
					EV::prod[n][t][i].set(GRB_DoubleAttr_Start, EV::val_prod[n][t][i]);
				}
			}
		}
		for (int b = 0; b < nBr; b++)
		{
			EV::Ze[b].set(GRB_DoubleAttr_Start, EV::val_Ze[b]);
		}
		Model.update();
	}

#pragma endregion

}




void NG_Module(GRBModel& Model, GRBLinExpr& exp_GVobj)
{

#pragma region Fetch Data
	vector<gnode> Gnodes = Params::Gnodes;
	vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	vector<eStore> Estorage = Params::Estorage;
	vector<exist_gSVL> Exist_SVL = Params::Exist_SVL;
	int nSVL = (int)Exist_SVL.size();
	vector<SVL> SVLs = Params::SVLs;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Enodes.size();
	int nPlt = (int)Plants.size();
	int nBr = (int)Branches.size();
	int neSt = (int)Estorage.size();
	int nPipe = (int)PipeLines.size();
	vector<int> Tg = Params::Tg;
	//vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	double pi = 3.14159;
	int nGnode = (int)Gnodes.size();
	double WACC = Params::WACC;
	int trans_unit_cost = Params::trans_unit_cost;
	int trans_line_lifespan = Params::trans_line_lifespan;
	int SVL_lifetime = Params::SVL_lifetime;
	double NG_price = Params::NG_price;
	double dfo_pric = Params::dfo_pric;
	double coal_price = Params::coal_price;
	double E_curt_cost = Params::E_curt_cost;
	double G_curt_cost = Params::G_curt_cost;
	double pipe_per_mile = Params::pipe_per_mile;
	int pipe_lifespan = Params::pipe_lifespan;
	map<int, vector<int>> Le = Params::Le;
	map<int, vector<int>> Lg = Params::Lg;
	vector<int> RepDaysCount = Params::RepDaysCount;
#pragma endregion

#pragma region NG network related costs
	//	GRBLinExpr exp_GVobj(env);
	GRBLinExpr exp_NG_import_cost(0);
	GRBLinExpr exp_invG(0);
	GRBLinExpr exp_pipe(0);
	GRBLinExpr exp_ng_curt(0);
	GRBLinExpr exp_rng_curt(0);
	GRBLinExpr exp_StrFOM(0);

	// cost of establishing new inter-network pipelines
	for (int i = 0; i < nPipe; i++)
	{
		double s1 = std::pow(1.0 / (1 + WACC), pipe_lifespan);
		double capco = WACC / (1 - s1);
		double cost1 = PipeLines[i].length * pipe_per_mile;

		exp_pipe += capco * cost1 * GV::Zg[i];
	}

	// fuel cost
	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			exp_NG_import_cost += RepDaysCount[tau] * GV::supply[k][tau] * NG_price;
		}
	}

	// storage investment cost
	for (int j = 0; j < nSVL; j++)
	{
		double s1 = (double)std::pow(1.0 / (1 + WACC), SVL_lifetime);
		double capco = WACC / (1 - s1);
		exp_invG += capco * (SVLs[0].Capex * GV::Xstr[j] + SVLs[1].Capex * GV::Xvpr[j]);
	}

	// storage FOM cost
	for (int j = 0; j < nSVL; j++)
	{
		exp_StrFOM += SVLs[0].FOM * (Exist_SVL[j].vap_cap + GV::Xvpr[j]);
		exp_StrFOM += SVLs[1].FOM * (Exist_SVL[j].store_cap + GV::Xstr[j]);
	}

	// Load shedding cost
	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			exp_ng_curt += RepDaysCount[tau] * G_curt_cost * GV::curtG[k][tau];
			exp_rng_curt += RepDaysCount[tau] * Params::RNG_price * GV::curtRNG[k][tau];
		}
	}


	exp_GVobj += exp_NG_import_cost + exp_invG + exp_pipe + exp_ng_curt + exp_rng_curt + exp_StrFOM;

	Model.addConstr(GV::strInv_cost == exp_invG);
	Model.addConstr(GV::pipe_cost == exp_pipe);
	Model.addConstr(GV::gShedd_cost == exp_ng_curt);
	Model.addConstr(GV::rngShedd_cost == exp_rng_curt);
	Model.addConstr(GV::gStrFOM_cost == exp_StrFOM);
	Model.addConstr(GV::NG_import_cost == exp_NG_import_cost);
	Model.addConstr(GV::NG_system_cost == exp_GVobj);
#pragma endregion

#pragma region NG Network Constraint
	// C1, C2: flow limit for NG
	for (int i = 0; i < nPipe; i++)
	{
		int is_exist = PipeLines[i].is_exist;
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			if (is_exist == 1)
			{
				Model.addConstr(GV::flowGG[i][tau] <= PipeLines[i].cap);
			}
			else
			{
				Model.addConstr(GV::flowGG[i][tau] <= PipeLines[i].cap * GV::Zg[i]);
			}
		}
	}

	//C3: flow balance, NG node
	double total_dem = 0;
	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			GRBLinExpr exp_gg_exp(0);
			GRBLinExpr exp_gg_imp(0);
			GRBLinExpr exp_ge_flow(0);
			GRBLinExpr exp_store(0);
			for (int l : Gnodes[k].Lexp) { exp_gg_exp += GV::flowGG[l][tau]; }
			for (int l : Gnodes[k].Limp) { exp_gg_imp += GV::flowGG[l][tau]; }

			for (int n : Gnodes[k].adjE) { exp_ge_flow += GV::flowGE[k][n][tau]; }
			for (int j : Gnodes[k].adjS)
			{
				exp_store += GV::flowVG[j][k][tau] - GV::flowGL[k][j][tau];
			}

			GRBLinExpr exp_curt(0);
			exp_curt = GV::curtG[k][tau] + GV::curtRNG[k][tau];
			double dem = Gnodes[k].demG[Tg[tau]];// +Gnodes[k].out_dem;
			total_dem += RepDaysCount[tau] * dem;
			//Model.addConstr(GV::supply[k][tau] + exp_gg_flow - exp_store + GV::curtG[k][tau] == Gnodes[k].demG[Tg[tau]] + Gnodes[k].out_dem);
			Model.addConstr(GV::supply[k][tau] - exp_gg_exp + exp_gg_imp - exp_ge_flow + exp_store + exp_curt == dem);
		}
	}


	// C3,C4: injection (supply) limit and curtailment limit
	GRBLinExpr exp_RNG(0);
	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			Model.addConstr(GV::supply[k][tau] >= Gnodes[k].injL);
			Model.addConstr(GV::supply[k][tau] <= Gnodes[k].injU);
			double dem = Gnodes[k].demG[Tg[tau]];// +Gnodes[k].out_dem;
			Model.addConstr(GV::curtG[k][tau] + GV::curtRNG[k][tau] <= dem);
			exp_RNG += RepDaysCount[tau] * GV::curtRNG[k][tau];
		}
	}
	Model.addConstr(exp_RNG <= Params::RNG_cap * Setting::PGC);

	//C5: storage balance
	for (int j = 0; j < nSVL; j++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			if (tau == 0)
			{
				//Model.addConstr(GV::Sstr[j][tau] == Exist_SVL[j].store_cap*0.5);
				// start with no stored gas
				Model.addConstr(GV::Sstr[j][tau] == Exist_SVL[j].store_cap * 0 + GV::Sliq[j][tau] - GV::Svpr[j][tau] / SVLs[1].eff_disCh);
				//Model.addConstr(GV::Sstr[j][tau] == GV::Sliq[j][tau] - GV::Svpr[j][tau] / SVLs[1].eff_disCh);
				continue;
			}
			Model.addConstr(GV::Sstr[j][tau] == (1 - SVLs[0].BOG) * GV::Sstr[j][tau - 1] + GV::Sliq[j][tau] - GV::Svpr[j][tau] / SVLs[1].eff_disCh);
		}
	}

	//C6: calculate Sliq
	for (int j = 0; j < nSVL; j++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			GRBLinExpr exp(0);
			// get the ng nodes adjacent to SVL j
			vector<int> NG_adj;
			for (int k = 0; k < nGnode; k++)
			{
				for (int j2 : Gnodes[k].adjS)
				{
					if (j2 == j)
					{
						NG_adj.push_back(k);
					}
				}
			}
			for (int k : NG_adj)
			{
				exp += GV::flowGL[k][j][tau];
			}
			Model.addConstr(GV::Sliq[j][tau] == exp);
		}
	}

	//C6: Sliq limit
	for (int j = 0; j < nSVL; j++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			Model.addConstr(GV::Sliq[j][tau] <= Exist_SVL[j].liq_cap);
		}
	}

	//C7: calculate Svpr
	for (int j = 0; j < nSVL; j++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			GRBLinExpr exp(0);
			// get the ng nodes adjacent to SVL j
			vector<int> NG_adj;
			for (int k = 0; k < nGnode; k++)
			{
				for (int j2 : Gnodes[k].adjS)
				{
					if (j2 == j)
					{
						NG_adj.push_back(k);
					}
				}
			}
			for (int k : NG_adj)
			{
				exp += GV::flowVG[j][k][tau];
			}
			Model.addConstr(GV::Svpr[j][tau] == exp);
		}
	}

	//C8: Svpr limit
	for (int j = 0; j < nSVL; j++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			Model.addConstr(GV::Svpr[j][tau] <= Exist_SVL[j].vap_cap + GV::Xvpr[j]);
		}
	}

	//C8: Sstr limit
	for (int j = 0; j < nSVL; j++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			Model.addConstr(GV::Sstr[j][tau] <= Exist_SVL[j].store_cap + GV::Xstr[j]);
		}
	}

#pragma endregion


#pragma region add the warm start solution
	if (Setting::warm_start_active && !Setting::heuristics1_active)
	{
		for (int i = 0; i < nPipe; i++)
		{
			GV::Zg[i].set(GRB_DoubleAttr_Start, GV::val_Zg[i]);
		}
		for (int j = 0; j < nSVL; j++)
		{
			GV::Xstr[j].set(GRB_DoubleAttr_Start, GV::val_Xstr[j]);
		}

		Model.update();
	}

#pragma endregion

}

void Coupling_Constraints(GRBModel& Model, GRBLinExpr& ex_xi, GRBLinExpr& ex_NG_emis, GRBLinExpr& ex_E_emis)
{

	CV::xi = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	CV::NG_emis = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	CV::E_emis = Model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);

#pragma region Fetch Data
	vector<gnode> Gnodes = Params::Gnodes;
	//vector<pipe> PipeLines = Params::PipeLines;
	vector<enode> Enodes = Params::Enodes;
	vector<plant> Plants = Params::Plants;
	//vector<eStore> Estorage = Params::Estorage;
	//vector<exist_gSVL> Exist_SVL = Params::Exist_SVL;
	//int nSVL = Exist_SVL.size();
	vector<SVL> SVLs = Params::SVLs;
	vector<branch> Branches = Params::Branches;
	int nEnode = (int)Enodes.size();
	int nPlt = (int)Plants.size();
	int nBr = (int)Branches.size();
	//int neSt = Estorage.size();
	//int nPipe = PipeLines.size();
	vector<int> Tg = Params::Tg;
	vector<int> Te = Params::Te;
	vector<int> time_weight = Params::time_weight;
	//double pi = 3.1415926535897;
	int nGnode = (int)Gnodes.size();
	map<int, vector<int>> Le = Params::Le;
	map<int, vector<int>> Lg = Params::Lg;
	vector<int> RepDaysCount = Params::RepDaysCount;
	double NG_emis_rate = Params::NG_emis_rate;
#pragma endregion


#pragma region Coupling 1
	//IloExpr ex_xi(env);
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
							if (Setting::DGSP_active)
							{
								//double s2 = std::max(0.001, EV::val_prod[n][t][i]);
								exp2 += time_weight[t] * Plants[i].heat_rate * EV::val_prod[n][t][i];
								flowge += time_weight[t] * Plants[i].heat_rate * EV::val_prod[n][t][i];
							}
							else
							{
								exp2 += time_weight[t] * Plants[i].heat_rate * EV::prod[n][t][i];
							}
						}
					}
				}
				Model.addConstr(RepDaysCount[tau] * GV::flowGE[k][n][tau] - exp2 <= 1e-2);
				Model.addConstr(exp2 - RepDaysCount[tau] * GV::flowGE[k][n][tau] <= 1e-2);
				ex_xi += RepDaysCount[tau] * GV::flowGE[k][n][tau];
			}
		}
	}
#pragma endregion

#pragma region Coupling 2
	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			ex_NG_emis += RepDaysCount[tau] * NG_emis_rate * (GV::supply[k][tau]);
		}
	}
	ex_NG_emis -= NG_emis_rate * CV::xi;

	for (int n = 0; n < nEnode; n++)
	{
		for (int t = 0; t < Te.size(); t++)
		{
			for (int i = 0; i < nPlt; i++)
			{
				if (Plants[i].type == "dfo" || Plants[i].type == "coal" || Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
				{
					ex_E_emis += time_weight[t] * Plants[i].emis_rate * Plants[i].heat_rate * EV::prod[n][t][i];
				}
			}
		}
	}
#pragma endregion

}
