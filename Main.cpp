// NGES.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include "ProblemData.h"
#include"Models_Funcs.h"
#pragma region Declaration of fields of "Setting" struct
bool Setting::print_E_vars = false;
bool Setting::print_NG_vars = false;
bool Setting::print_results_header = true;
bool Setting::relax_int_vars = false;
bool Setting::relax_UC_vars = false;
bool Setting::is_xi_given = false;
double Setting::xi_val;
bool Setting::xi_UB_obj = false;
bool Setting::xi_LB_obj = false;
double Setting::cplex_gap;
double Setting::CPU_limit;
int Setting::Num_rep_days;
double Setting::Emis_lim;
double Setting::RPS;
bool Setting::DGSP_active = false;
bool Setting::DESP_active = false;
bool Setting::Approach_1_active = true;
bool Setting::Approach_2_active = false;
double Setting::RNG_cap;
int Setting::Case;
//bool Setting::MP_init_heuristic = false;
bool Setting::warm_start_active = false;
bool Setting::print_all_vars = false;
bool Setting::fix_some_E_NG_vars = false;
double Setting::PGC;
double Setting::PE;
bool Setting::use_benders = false;
bool Setting::UC_active = false;
bool Setting::multi_cut_active = false;
#pragma endregion

int main(int argc, char* argv[])
{
	auto start = chrono::high_resolution_clock::now();
	GRBEnv* env; env = new GRBEnv();
	//double total_ng_inj_per_year = 4.9E+09; //MMBtu
	//double total_possible_emis_per_year = total_ng_inj_per_year * 0.05831;//=2.86e8
	double emission_power_plants = 24e6; // 24 tons = 24e9 kg
	double total_ng_yearly_demand = 5.49E+08;// MMBtu in 2018
	double total_yearly_gen_demand = emission_power_plants / 0.053; //4.52e8 MMBtu
	double PGC = total_ng_yearly_demand + total_yearly_gen_demand;//1.00E+09 MMBTu possible gas consumption
	double Poss_Emis = PGC * 0.053; //5.3e7 tons  (previsouly 8.2e7 kg)
	// 0.053 tons of co2 per MMBTus of NG
	Setting::PGC = PGC;
	Setting::PE = Poss_Emis;
	if (argc > 1)
	{
		Setting::Num_rep_days = atoi(argv[1]);
		Setting::Approach_1_active = atoi(argv[2]);
		Setting::Approach_2_active = atoi(argv[3]);
		Setting::Case = atoi(argv[4]);
		Setting::is_xi_given = atoi(argv[5]);
		Setting::xi_val = atof(argv[6]);
		Setting::Emis_lim = atof(argv[7]);
		Setting::RPS = atof(argv[8]);
		Setting::RNG_cap = atof(argv[9]);
		Setting::cplex_gap = atof(argv[10]);  
		Setting::CPU_limit = atoi(argv[11]);   // seconds
		Setting::UC_active = atoi(argv[12]);
		Setting::relax_UC_vars = atoi(argv[13]);
		Setting::use_benders = atoi(argv[14]);
		Setting::multi_cut_active = atoi(argv[15]);
	}
	else
	{
		Setting::Num_rep_days = 30;   // 2,5, 7,10 14,15,20 30, 52,104, 365
		Setting::Approach_1_active = true; // approach 1: integrated, 2: decoupled 
		Setting::Approach_2_active = false; // default = false
		Setting::Case = 3; //1: indep. networks, 2: only E emission, 3:joint planning
		Setting::is_xi_given = false;
		Setting::xi_val = 0.0;//0.01,0.05, 0.1,0.15,0.2,;
		Setting::Emis_lim = 0.2;    // xPE= tons  (for case 2: 9%PE~20% of EE (elec emission),  
		Setting::RPS = 0.6;		    // out of 1 (=100%) Renewable Portfolio Share
		Setting::RNG_cap = 0.4; //0.2,0.3,0.4,
		Setting::cplex_gap = 0.01;  // 
		Setting::CPU_limit = 10800;   // seconds		
		Setting::UC_active = true; // only applies to the full problem
		Setting::relax_UC_vars = false;
		Setting::relax_int_vars = false; // int vars (but not binary vars) in electricity network
		Setting::use_benders = true;
		Setting::multi_cut_active = false;
	}

#pragma region Problem Setting
	Setting::print_results_header = false;
	//Setting::xi_LB_obj = false; // (default = false) 
	//Setting::xi_UB_obj = false;  // (default = false) 

	bool only_feas_sol = false;
	Setting::print_all_vars = true;

#pragma endregion

#pragma region  Other parameters   
	Params::Num_Rep_Days = Setting::Num_rep_days;
	double RNG_price = 20; // $$ per MMBtu
	double WACC = 0.071;// Weighted average cost of capital to calculate CAPEX coefficient from ATB2021
	int trans_unit_cost = 3500; // dollars per MW per mile of trans. line (ReEDS 2019)
	int trans_line_lifespan = 30; // years
	int decom_lifetime = 2035 - 2016;
	int battery_lifetime = 15; // 
	double NG_price = 5.45;//per MMBTu, approximated from NG price in eia.gov
	double dfo_pric = (1e6 / 1.37e5) * 3.5;//https://www.eia.gov/energyexplained/units-and-calculators/ and https://www.eia.gov/petroleum/gasdiesel/
	double coal_price = 92 / 19.26; //https://www.eia.gov/coal/ and https://www.eia.gov/tools/faqs/faq.php?id=72&t=2#:~:text=In%202020%2C%20the%20annual%20average,million%20Btu%20per%20short%20ton.
	double Nuclear_price = 0.72; // per MMBtu from 2045 ATB 2021
	double E_curt_cost = 10e3; // $ per MWh;
	double G_curt_cost = 1e3; // & per MMBtu
	double pipe_per_mile = 7e+5;//https://www.gem.wiki/Oil_and_Gas_Pipeline_Construction_Costs
	int SVL_lifetime = 30; //https://www.hydrogen.energy.gov/pdfs/19001_hydrogen_liquefaction_costs.pdf
	int pipe_lifespan = 30; // years, https://www.popsci.com/story/environment/oil-gas-pipelines-property/#:~:text=There%20are%20some%203%20million,%2C%20power%20plants%2C%20and%20homes.&text=Those%20pipelines%20have%20an%20average%20lifespan%20of%2050%20years.
	//double Ng_demand_growth_by_2050 = 0.5; // 50% https://www.eia.gov/todayinenergy/detail.php?id=42342
	double NG_emis_rate = 0.05331;  // tons of CO2 per MMBtu
#pragma endregion

#pragma region Read Data
	Read_Data();
#pragma endregion

#pragma region Set Params
	Params::WACC = WACC;
	Params::trans_unit_cost = trans_unit_cost;
	Params::trans_line_lifespan = trans_line_lifespan;
	Params::NG_price = NG_price;
	Params::dfo_pric = dfo_pric;
	Params::coal_price = coal_price;
	Params::nuclear_price = Nuclear_price;
	Params::E_curt_cost = E_curt_cost;
	Params::G_curt_cost = G_curt_cost;
	Params::pipe_per_mile = pipe_per_mile;
	Params::pipe_lifespan = pipe_lifespan;
	Params::SVL_lifetime = SVL_lifetime;
	Params::battery_lifetime = battery_lifetime;
	Params::RNG_cap = Setting::RNG_cap;
	Params::RNG_price = RNG_price;
	//Params::RPS = RPS;
	Params::NG_emis_rate = NG_emis_rate;

#pragma endregion
	int nBr = Params::Branches.size();

	int** Xs = new int* [Params::Enodes.size()];
	double*** Ps = new double** [Params::Enodes.size()];
	int** XestS = new int* [Params::Enodes.size()];
	int** XdecS = new int* [Params::Enodes.size()];


	double opt = 0; double UB1 = 0; double UB2 = 0;
	double MidSol = 0;
	double xiLB1 = 0; double xiUB1 = 0;

	double elec_LB = 0;
	double elec_UB = 0;
	double ng_obj;
	double feas_gap;



	if (Setting::use_benders)
	{
		std::cout << "using Benders..." << endl;
		double total_cost = Benders_Decomposition(env);
		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds

		Print_Results(Elapsed, total_cost);
		return 0;
	}


	if (only_feas_sol)
	{
		double UB = feas_sol(elec_LB, elec_UB, ng_obj, feas_gap);
		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		//Print_Results(Elapsed, EV::val_e_system_cost + GV::val_NG_system_cost);
		std::cout << "\n\n \t\t Elapsed time: " << Elapsed;
		std::cout << "\t\t System LB: " << elec_LB + ng_obj;
		std::cout << "\t\t System UB: " << elec_UB + ng_obj << endl;
		std::cout << "\n\n";
		return 0;
	}


	if (Setting::Case == 1)
	{
		Electricy_Network_Model(env);
		NG_Network_Model(env);
		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		Print_Results(Elapsed, 0);
		return 0;
	}


	if (Setting::Approach_1_active)
	{
		std::cout << "case 3 by approach 1" << endl;
		bool ap2 = Setting::Approach_2_active;
		Setting::Approach_2_active = false;
		double total_cost = Integrated_Model(env,Setting::cplex_gap);
		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		Print_Results(Elapsed, total_cost);
		if (Setting::xi_LB_obj) { xiLB1 = Integrated_Model(env, Setting::cplex_gap); }
		if (Setting::xi_UB_obj) { xiUB1 = Integrated_Model(env, Setting::cplex_gap); }
		start = chrono::high_resolution_clock::now();
		Setting::Approach_2_active = ap2;
	}





	double desp_obj = 0;
	double dgsp_obj = 0; double xiLB2 = 0; double xiUB2 = 0;
	if (Setting::Approach_2_active)
	{
		Setting::Approach_1_active = false;
		desp_obj = DESP(env);
		dgsp_obj = DGSP(env);
		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		Print_Results(Elapsed, dgsp_obj);
	}


}
