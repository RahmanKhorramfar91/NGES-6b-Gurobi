#pragma once
#include <iostream>
#include <vector> 
#include<random>
#include <string>
#include <fstream> // to read and write from/on files
#include <cstdio>  // for puts and printf
#include <cstdlib> // for string converstion, rand srand, dynamic memory management (C++ header)
#include <iterator>
#include<sstream>
#include<chrono> // for high resolution timing (elapsed time in microseconds and so on)
#include<array>
#include<list>
#include<cmath> // power operator and others
#include<numeric>
#include<map>  // To define dictionary data types

using namespace std;
void Read_rep_days(string name, vector<int>& Rep, vector<int>& RepCount);
struct enode
{
	int num;
	//int fips;
	//int fips_order;
	vector<int> adj_buses;
	//vector<double> adj_dist;
	vector<string> Init_plt_types;
	vector<int> Init_plt_ind;
	vector<double*> Init_plt_prod_lim; // [Pmin,Pmax]
	vector<int> Init_plt_count;
	vector<double> demand;
	enode(int n)
	{
		this->num = n;
	}

	static vector<enode> read_bus_data(string FileName);
	static void read_adj_data(string FileName, vector<enode>& Enodes);
	static void read_exist_plt_data(string FileName, vector<enode>& Enodes);
	static void read_demand_data(string FileName, vector<enode>& Enodes);
};

//struct Fips
//{
//	int fips_num;
//	int bus_count;
//
//	vector<double> demE;
//	Fips(int f, int bc)
//	{
//		this->fips_num = f;
//		this->bus_count = bc;
//	}
//
//	static vector<Fips> read_fips_data(string FileName);
//	static void read_demand_Data(string FileName, vector<Fips>& all_FIPS);
//};

struct plant
{
	string type;
	int num;  //0:'ng',1:'dfo',2:'solar',3:'wind',4:'wind_offshore', 5:'hydro',6:'coal',7:'nuclear'
	int is_exis; // 1=the plant existing type (PowerSim Data), 0= the plant is new type (ATB)
	//int Umax; // maximum number of plants in a Enode
	double capex; // dollar per plant
	int fix_cost;//dollar per count of type i per year
	int var_cost;//dollar per MWh
	double emis_rate;//ton/MMBTU
	double heat_rate;//MMBTU/MWh
	int lifetime; // year
	int decom_cost;
	double Pmax;  // MW
	double Pmin;
	double rampU;
	//double rampD;
	//int emis_cost; //$/ton
	// revisit these paramters laters


	vector<vector<double>> zonal_profile;
	vector<double> Reg_coeffs_per_state;

	plant(string t, int n, int ise, int f, int v, double emi, double hr, int lt, int dec,
		double pmax, double pmin, double ru, double cap)
	{
		this->type = t;
		this->num = n;
		this->is_exis = ise;
		this->capex = cap;
		this->fix_cost = f;
		this->var_cost = v;
		this->emis_rate = emi;
		this->heat_rate = hr;
		this->lifetime = lt;
		this->decom_cost = dec;
		this->Pmax = pmax;
		this->Pmin = pmin;
		this->rampU = ru;
		//this->rampD = rd;
		//this->emis_cost = emic;
		//this->Umax = max_num;
	}

	static vector<plant> read_new_plant_data(string FileName);
	static vector<plant> read_regional_coeffs(string FileName, vector<plant>& Plants);
	static void read_VRE_profile(string FileName1, string FileName2, string FileName3, vector<plant>& Plants);
};

struct branch
{
	int from_bus;
	int to_bus;
	double suscep;
	double maxFlow;
	double length;
	int is_exist;
	branch(int f, int t, double l, int ie, double m, double s)
	{
		this->from_bus = f;
		this->to_bus = t;
		this->length = l;
		this->is_exist = ie;
		this->maxFlow = m;
		this->suscep = s;
	}
	static vector<branch> read_branch_data(int nBus, string FN,
		map<int, vector<int>>& Le);

};

struct eStore
{
	int energy_cost;// &/MWh
	int power_cost; // &/MW
	double eff_ch;  // charge efficiency
	double eff_disCh; // discharge efficiency
	double eFOM; // energy related fixed operating and maintanence cost
	double pFOM; // power related fixed operating and maintanence cost



	eStore(int en, int pow, double ch, double dis, double efom, double pfom)
	{
		this->energy_cost = en;
		this->power_cost = pow;
		this->eff_ch = ch;
		this->eff_disCh = dis;
		this->eFOM = efom;
		this->pFOM = pfom;
	}
	static vector<eStore> read_elec_storage_data(string FileName);

};

struct gnode
{
	int num;
	int fips;
	vector<double> demG;
	int out_dem;
	int injU;
	int injL = 0;
	vector<int> Lexp;
	vector<int> Limp;
	vector<int> adjE;
	vector<int> adjS;
	gnode(int n, int f, int od, int inj, vector<int> as)
	{
		this->num = n;
		this->fips = f;
		this->out_dem = od;
		this->injU = inj;
		this->adjS = as;
	}

	static vector<gnode> read_gnode_data(string FileName);
	static void read_Lexp_data(string FileName, vector<gnode>& Gnodes);
	static void read_Limp_data(string FileName, vector<gnode>& Gnodes);
	static void read_adjE_data(string FileName, vector<gnode>& Gnodes);
	static void read_ng_demand_data(string FileName, vector<gnode>& Gnodes);
};

struct pipe
{
	int from_node;
	int to_node;
	int is_exist;
	double length;
	double cap;

	pipe(int f, int t, int exi, double le, double ca)
	{
		this->from_node = f;
		this->to_node = t;
		this->is_exist = exi;
		this->length = le;
		this->cap = ca;
	}
	static vector<pipe> read_pipe_data(string FileName, map<int, vector<int>>& Lg);
};

struct exist_gSVL
{
	int num;
	int store_cap; // max storage capacity
	int vap_cap; // max vaporization capacity
	int liq_cap; // max liquefaction capacity

	exist_gSVL(int n, int sc, int vc, int lc)
	{
		this->num = n;
		this->store_cap = sc;
		this->vap_cap = vc;
		this->liq_cap = lc;
	}
	static vector<exist_gSVL> read_exist_SVL_data(string FL);
};

struct SVL
{
	// SVL 0 : Storage facility
	// SVL 1: Vaporization facility
	// SVL 2: Liquefaction facility
	int Capex;
	int FOM;
	double eff_ch;  // charge efficiency
	double eff_disCh; // discharge efficiency
	double BOG = 6e-6; // Boil-off gas, Dharik: The boiled off gas is sent to the pipeline network. You can convert this to an hourly rate and model this as the minimum discharge of the tank in each hour
	SVL(int cp, int fm, double ch, double dis)
	{
		this->Capex = cp;
		this->FOM = fm;
		this->eff_ch = ch;
		this->eff_disCh = dis;
	}

	static vector<SVL> read_SVL_data(string FL);
};

struct Params
{
	// General
	static int Num_Rep_Days;
	static double WACC;
	static int trans_unit_cost;
	static int trans_line_lifespan;
	static double NG_price;
	static double RNG_price;
	static double RNG_cap;

	static double dfo_pric;
	static double coal_price;
	static double nuclear_price;
	static double E_curt_cost;
	static double G_curt_cost;
	static double pipe_per_mile;
	static int pipe_lifespan;
	static int battery_lifetime;
	static int SVL_lifetime;
	static vector<int> RepDaysCount;
	static double NG_emis_rate;


	// natural gas data
	static vector<int> Tg;
	static vector<gnode> Gnodes;
	static vector<pipe> PipeLines;
	static map<int, vector<int>> Lg;
	static vector<exist_gSVL> Exist_SVL;
	static vector<SVL> SVLs;


	// electricity data
	static map<int, vector<int>> Le;
	static vector<int> Te;
	static vector<int> time_weight;
	static vector<enode> Enodes;
	static vector<plant> Plants;
	static vector<branch> Branches;
	static vector<eStore> Estorage;
};


void Read_Data();


