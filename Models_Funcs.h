#pragma once
#include"ProblemData.h"
#include "gurobi_c++.h"

void Electricy_Network_Model();
double feas_sol(double& elec_LB, double& elec_UB, double& ng_obj, double& feas_gap);
void NG_Network_Model();




double Integrated_Model();

double DGSP();
double DESP();


void Read_rep_days(string name, vector<int>& Rep, vector<int>& RepCount);

void Populate_EV_SP(GRBModel& Model);
void Elec_Module(GRBModel& Model, GRBLinExpr& exp_Eobj);
void Populate_GV(GRBModel& Model);
void NG_Module(GRBModel& Model, GRBLinExpr& exp_GVobj);

void Populate_GV(GRBModel& Model);



void Coupling_Constraints(GRBModel& Model, GRBLinExpr& ex_xi, GRBLinExpr& ex_NG_emis, GRBLinExpr& ex_E_emis);

void Get_EV_vals(GRBModel Model);
void Get_GV_vals(GRBModel Model);

void Print_Results(double Elapsed_time, double status);


struct Setting
{
	static bool print_E_vars;
	static bool print_NG_vars;
	static bool relax_int_vars;
	static bool print_results_header;
	static bool heuristics1_active;
	static bool warm_start_active;

	static bool is_xi_given;
	static double xi_val;
	static bool xi_UB_obj;
	static bool xi_LB_obj;
	static double cplex_gap;
	static double CPU_limit;

	static int Num_rep_days;
	static double Emis_lim;
	static double RPS;
	static bool Approach_1_active;
	static bool Approach_2_active;
	static int Case;

	static bool DGSP_active;
	static bool DESP_active;
	static double RNG_cap;
	static bool print_all_vars;

	static double PGC;
	static double PE;
};

struct EV
{
	static 	GRBVar** Xest; // integer (continues) for plants
	static	GRBVar** Xdec; // integer (continues) for plants
	static	GRBVar** YeCD; // continuous: charge/discharge capacity
	static	GRBVar** YeLev; // continuous: charge/discharge level
	static GRBVar** YeStr; //
	static	GRBVar* Ze;// binary
	static	GRBVar** theta; // continuous phase angle
	static	GRBVar** curtE; // continuous curtailment variable
	static	GRBVar*** prod;// continuous
	static	GRBVar*** eSch;// power charge to storage 
	static	GRBVar*** eSdis;// power discharge to storage
	static	GRBVar*** eSlev;// power level at storage
	static	GRBVar** Xop; // integer (continues)
	static	GRBVar** flowE; // flowE subscripts are "ntm" here
	static	GRBVar est_cost; // cost of establishments
	static	GRBVar decom_cost;
	static	GRBVar fixed_cost;
	static	GRBVar var_cost;
	static	GRBVar thermal_fuel_cost;
	static	GRBVar shedding_cost;
	static	GRBVar elec_storage_cost;
	static	GRBVar Emit_var;
	static	GRBVar dfo_coal_emis_cost;
	static GRBVar e_system_cost;



	static double*** val_prod;
	static double*** val_eSch;
	static double*** val_eSdis;
	static	double val_est_cost;
	static	double val_decom_cost;
	static	double val_fixed_cost;
	static	double val_var_cost;
	static	double val_thermal_fuel_cost;
	static	double val_shedding_cost;
	static	double val_elec_storage_cost;
	static	double val_Emit_var;
	static double val_dfo_coal_emis_cost;
	static  double val_e_system_cost;
	static double val_num_storage; // number of storage facilities established
	static double val_storage_lev; // value of toal storage capacity (MWh)
	static double val_storage_cap; // value of toal storage capacity (MW)
	static double* val_num_est;// number of established plants
	static double* val_num_decom; // number of plants decommissioned
	static double val_total_curt; // total load shedding
	static double val_num_est_trans; // number of new transmission lines established
	static double val_total_flow; // total flow in the E network
	static double* val_total_prod; // total production	
	static double** val_Xop;
	static double** val_Xest;
	static double** val_Xdec;
	static double* val_Ze;
	static double** val_YeStr;
	static double MIP_gap;
	static double** val_curtE;
	static double*** val_eSlev;



};

struct GV
{

	static GRBVar* Xstr;
	static GRBVar* Xvpr;
	static GRBVar** Sstr;
	static GRBVar** Svpr;
	static GRBVar** Sliq;
	static GRBVar** supply;
	static GRBVar** curtG;
	static GRBVar** curtRNG;
	static GRBVar** flowGG;
	static GRBVar*** flowGE;
	static GRBVar*** flowGL;
	static GRBVar*** flowVG;
	static GRBVar* Zg;

	static GRBVar strInv_cost;
	static GRBVar pipe_cost;
	static GRBVar gShedd_cost;
	static GRBVar rngShedd_cost;
	static GRBVar gStrFOM_cost;
	static GRBVar NG_import_cost;
	static GRBVar NG_system_cost;

	static double*** val_flowGE;
	static double val_strInv_cost;
	static double val_pipe_cost;
	static double val_ngShedd_cost;
	static double val_rngShedd_cost;
	static double val_gStrFOM_cost;
	static double val_NG_import_cost;
	static double val_NG_system_cost;

	static double* val_total_nodal_supply;
	static double** val_supply;
	static double val_ng_curt;
	static double val_rng_curt;
	static int val_num_est_pipe;
	static double val_total_flowGG;
	static double val_total_flowGE;
	static double val_total_flowGL;
	static double val_total_flowVG;
	static double* val_storage;
	static double* val_vapor;
	static double* val_Xstr;
	static double* val_Zg;
	static double MIP_gap;
};

struct CV
{
	static GRBVar xi; // flow from NG to E network (NG consumed by NG-fired plants)
	static GRBVar NG_emis; // emission from NG network
	static GRBVar E_emis; // emission from Electricity network. 


	static double used_emis_cap; // eta value from the DGSP of the Approach 2

	static double val_xi; // f
	static double val_NG_emis;
	static double val_E_emis;

};


void Subproblem();
void Primal_subproblem();
void Elec_Module_Primal_SP(GRBModel& Model, GRBLinExpr& exp_Eobj);
struct MP
{
	static GRBVar Psi;
};

struct SP
{
	static double SP_Primal_obj;
	static double SP_Dual_obj;
	// Dual variables for power balance constraint to get the locational price of natural gas
// or to get dual information to add cuts
	static GRBConstr* d_alpha;
	static GRBConstr*** d_beta;
	static GRBConstr*** d_gamma1;
	static GRBConstr*** d_gamma2;
	static GRBConstr** d_delta11;
	static GRBConstr** d_delta12;
	static GRBConstr** d_delta21;
	static GRBConstr** d_delta22;
	static GRBConstr** d_theta;// = new GRBConstr * [Params::Plants.size()];
	static GRBConstr** d_zeta11;
	static GRBConstr** d_zeta12;
	static GRBConstr** d_zeta21;
	static GRBConstr** d_zeta22;
	static GRBConstr** d_eta1;
	static GRBConstr** d_eta2;
	static GRBConstr*** d_pi;
	static GRBConstr** d_phi;
	static GRBConstr d_omega;
	static GRBConstr*** d_rho;
	static GRBConstr d_tau;

	static double* dual_val_alpha;
	static double*** dual_val_beta;
	static double*** dual_val_gamma1;
	static double*** dual_val_gamma2;
	static double** dual_val_delta11;
	static double** dual_val_delta12;
	static double** dual_val_delta21;
	static double** dual_val_delta22;
	static double** dual_val_theta; // power balance equation dual variables
	static double** dual_val_zeta11;
	static double** dual_val_zeta12;
	static double** dual_val_zeta21;
	static double** dual_val_zeta22;
	static double** dual_val_eta1;
	static double** dual_val_eta2;
	static double*** dual_val_pi;
	static double** dual_val_phi;
	static double dual_val_omega;
	static double*** dual_val_rho;
	static double dual_val_tau;



	// from coupling constraints
	static GRBVar*** rho;
	static GRBVar tau;

	static GRBVar* alpha;
	static GRBVar*** beta;
	static GRBVar*** gamma1;
	static GRBVar*** gamma2;
	static GRBVar** delta11;	static GRBVar** delta12;
	static GRBVar** delta21; 	static GRBVar** delta22;

	static GRBVar** theta;

	static GRBVar** zeta11;
	static GRBVar** zeta12;
	static GRBVar** zeta21;
	static GRBVar** zeta22;
	static GRBVar** eta1;
	static GRBVar** eta2;
	static GRBVar* eta3;

	static GRBVar*** pi;
	static GRBVar** phi;
	static GRBVar omega;
};