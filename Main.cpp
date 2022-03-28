// NGES.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include "ProblemData.h"
#include"Models_Funcs.h"
#pragma region Declaration of fields of "Setting" struct
bool Setting::print_E_vars;
bool Setting::print_NG_vars;
bool Setting::print_results_header;
bool Setting::relax_int_vars;
bool Setting::is_xi_given;
double Setting::xi_val;
bool Setting::xi_UB_obj;
bool Setting::xi_LB_obj;
double Setting::cplex_gap;
double Setting::CPU_limit;
int Setting::Num_rep_days;
double Setting::Emis_lim;
double Setting::RPS;
bool Setting::DGSP_active;
bool Setting::DESP_active;
bool Setting::Approach_1_active;
bool Setting::Approach_2_active;
double Setting::RNG_cap;
int Setting::Case;
bool Setting::heuristics1_active;
bool Setting::warm_start_active;
#pragma endregion

int main(int argc, char* argv[])
{
	auto start = chrono::high_resolution_clock::now();
	double total_ng_inj_per_year = 4.9E+09; //MMBtu
	double total_possible_emis_per_year = total_ng_inj_per_year * 0.05831;//=2.86e8
	double total_ng_yearly_demand = 7.02E+08;// MMBtu
	double total_yearly_gen_demand = 4.4e8; //MMBtu
	double PGC = total_ng_yearly_demand + total_yearly_gen_demand;//1.142E+09 possible gas consumption
	double Poss_Emis = PGC * 0.053; //8.2e7
	//440834343

	if (argc > 1)
	{
		Setting::Num_rep_days = atoi(argv[1]);
		Setting::Approach_1_active = atoi(argv[2]);
		Setting::Approach_2_active = atoi(argv[3]);
		Setting::Case = atoi(argv[4]);
		Setting::is_xi_given = atoi(argv[5]);
		Setting::xi_val = atof(argv[6]) * PGC;
		Setting::Emis_lim = atof(argv[7]) * Poss_Emis;
		Setting::RPS = atof(argv[8]);
		Setting::RNG_cap = atof(argv[9]) * PGC;
		Setting::cplex_gap = atof(argv[10]);  // 1%
		Setting::CPU_limit = atoi(argv[11]);   // seconds
	}
	else
	{
		Setting::Num_rep_days = 7;   // 2, 7, 14, 52, 365
		Setting::Approach_1_active = false; // approach 1: integrated, 2: decoupled 
		Setting::Approach_2_active = false; // default = false
		Setting::Case = 1; //1: indep. networks, 2: only E emission, 3:joint planning
		Setting::is_xi_given = false;
		Setting::xi_val = 0.1 * PGC;//0.01,0.05, 0.1,0.15,0.2,;
		Setting::Emis_lim = 0.55 * Poss_Emis;    // tons
		Setting::RPS = 0.1;		    // out of 1 (=100%) Renewable Portfolio Share
		Setting::RNG_cap = 0.2 * PGC; //0.05,0.1,0.2,
		Setting::cplex_gap = 0.01;  // 2%
		Setting::CPU_limit = 3600;   // seconds
	}

#pragma region Problem Setting
	Setting::relax_int_vars = false; // int vars in electricity network
	Setting::print_results_header = false;
	Setting::xi_LB_obj = false; // (default = false) 
	Setting::xi_UB_obj = false;  // (default = false) 
	//Setting::heuristics1_active = true;
	Setting::warm_start_active = false;
#pragma endregion

#pragma region  Other parameters   
	Params::Num_Rep_Days = Setting::Num_rep_days;
	double RNG_price = 20; // $$ per MMBtu
	double WACC = 0.05;// Weighted average cost of capital to calculate CAPEX coefficient from ATB2021
	int trans_unit_cost = 3500; // dollars per MW per mile of trans. line (ReEDS 2019)
	int trans_line_lifespan = 40; // years
	int decom_lifetime = 2035 - 2016;
	int battery_lifetime = 15; // 
	double NG_price = 4;//per MMBTu, approximated from NG price in eia.gov
	double dfo_pric = (1e6 / 1.37e5) * 3.5;//https://www.eia.gov/energyexplained/units-and-calculators/ and https://www.eia.gov/petroleum/gasdiesel/
	double coal_price = 92 / 19.26; //https://www.eia.gov/coal/ and https://www.eia.gov/tools/faqs/faq.php?id=72&t=2#:~:text=In%202020%2C%20the%20annual%20average,million%20Btu%20per%20short%20ton.
	double Nuclear_price = 0.69; // per MMBtu
	double E_curt_cost = 3e4; // $ per MWh;
	double G_curt_cost = 3e3; // & per MMBtu
	double pipe_per_mile = 7e+5;//https://www.gem.wiki/Oil_and_Gas_Pipeline_Construction_Costs
	int SVL_lifetime = 40; //https://www.hydrogen.energy.gov/pdfs/19001_hydrogen_liquefaction_costs.pdf
	int pipe_lifespan = 50; // years, https://www.popsci.com/story/environment/oil-gas-pipelines-property/#:~:text=There%20are%20some%203%20million,%2C%20power%20plants%2C%20and%20homes.&text=Those%20pipelines%20have%20an%20average%20lifespan%20of%2050%20years.
	double Ng_demand_growth_by_2050 = 0.5; // 50% https://www.eia.gov/todayinenergy/detail.php?id=42342
	double NG_emis_rate = 0.05831;  // tons of CO2 per MMBtu
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


	int** Xs = new int* [Params::Enodes.size()];
	double*** Ps = new double** [Params::Enodes.size()];
	int** XestS = new int* [Params::Enodes.size()];
	int** XdecS = new int* [Params::Enodes.size()];


	double opt = 0; double UB1 = 0; double UB2 = 0;
	double MidSol = 0;
	double xiLB1 = 0; double xiUB1 = 0;

	if (Setting::Case == 1)
	{
		Electricy_Network_Model();
		NG_Network_Model();
		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		Print_Results(Elapsed, EV::val_e_system_cost + GV::val_NG_system_cost);
	}


	if (Setting::Approach_1_active)
	{
		bool ap2 = Setting::Approach_2_active;
		Setting::Approach_2_active = false;
		double total_cost = Integrated_Model();
		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		Print_Results(Elapsed, total_cost);
		if (Setting::xi_LB_obj) { xiLB1 = Integrated_Model(); }
		if (Setting::xi_UB_obj) { xiUB1 = Integrated_Model(); }
		start = chrono::high_resolution_clock::now();
		Setting::Approach_2_active = ap2;
	}
	double desp_obj = 0;
	double dgsp_obj = 0; double xiLB2 = 0; double xiUB2 = 0;
	if (Setting::Approach_2_active)
	{
		Setting::Approach_1_active = false;
		desp_obj = DESP();
		dgsp_obj = DGSP();
		auto end = chrono::high_resolution_clock::now();
		double Elapsed = (double)chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000; // seconds
		Print_Results(Elapsed, dgsp_obj);
	}


}
