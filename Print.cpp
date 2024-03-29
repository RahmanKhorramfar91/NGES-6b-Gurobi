#include"Models_Funcs.h"

void Get_EV_vals(GRBModel Model)
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
#pragma endregion

#pragma region Get dual variables
	//EV::dual_val_theta = new double* [nEnode];
	//EV::dual_val_phi = new double* [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		//EV::dual_val_theta[n] = new double[Te.size()]();
		//EV::dual_val_phi[n] = new double[Te.size()]();
		double ave_price = 0;
		for (int t = 0; t < Te.size(); t++)
		{
			if (!Setting::relax_int_vars)
			{
				break;
			}
			//EV::dual_val_theta[n][t] = EV::d_theta[n][t].get(GRB_DoubleAttr_Pi);
			//EV::dual_val_phi[n][t] = EV::d_phi[n][t].get(GRB_DoubleAttr_Pi);
			//ave_price += SP::dual_val_theta1[n][t];
			/*if (EV::val_PB[n][t]>0)
			{
				cout << "PB[" << n << "][" << t << "] = " << EV::val_PB[n][t] << endl;
			}*/
		}
		//ave_price = ave_price / Te.size();
		//cout << "NG price at node \t" << n << ": \t" << ave_price << endl;
	}
#pragma endregion


#pragma region get Electricity Network variables
	//int periods2print = 20;
	/*ofstream fid;
	if (Setting::DESP_active)
	{
		fid.open("DESP.txt");
	}
	else
	{
		fid.open("DVe.txt");
	}*/
	EV::val_est_cost = EV::est_cost.get(GRB_DoubleAttr_X);
	EV::val_est_trans_cost = EV::est_trans_cost.get(GRB_DoubleAttr_X);
	EV::val_decom_cost = EV::decom_cost.get(GRB_DoubleAttr_X);
	EV::val_fixed_cost = EV::fixed_cost.get(GRB_DoubleAttr_X);
	EV::val_startup_cost = EV::startup_cost.get(GRB_DoubleAttr_X);
	EV::val_var_cost = EV::var_cost.get(GRB_DoubleAttr_X);
	EV::val_thermal_fuel_cost = EV::thermal_fuel_cost.get(GRB_DoubleAttr_X);
	EV::val_shedding_cost = EV::shedding_cost.get(GRB_DoubleAttr_X);
	EV::val_elec_storage_cost = EV::elec_storage_cost.get(GRB_DoubleAttr_X);

	EV::val_dfo_coal_emis_cost = EV::dfo_coal_emis_cost.get(GRB_DoubleAttr_X);
	EV::val_e_system_cost = EV::e_system_cost.get(GRB_DoubleAttr_X);

	if (Setting::Approach_1_active || Setting::Approach_2_active)
	{
		EV::val_Emit_var = CV::E_emis.get(GRB_DoubleAttr_X);
		CV::val_E_emis = CV::E_emis.get(GRB_DoubleAttr_X);
		CV::val_xi = CV::xi.get(GRB_DoubleAttr_X);
	}



	//fid << "Elapsed time: " << Elapsed_time << endl;
	//fid << "\t Total cost for both networks:" << obj << endl;
	//fid << "\t Electricity Network Obj Value:" << EV::e_system_cost) << endl;
	//fid << "\t Gap: " << gap << " Status:" << cplex.getStatus() << endl;
	//fid << "\n \t Establishment Cost: " << EV::est_cost);
	//fid << "\n \t Decommissioning Cost: " << EV::decom_cost);
	//fid << "\n \t Fixed Cost: " << EV::fixed_cost);
	//fid << "\n \t Variable Cost: " << EV::var_cost);
	////fid << "\n \t Emission Cost: " << emis_cost);
	//fid << "\n \t (dfo, coal, and nuclear) Fuel Cost: " << EV::thermal_fuel_cost);
	//fid << "\n \t Load Shedding Cost: " << EV::shedding_cost);
	//fid << "\n \t Storage Cost: " << EV::elec_storage_cost) << "\n";
	//fid << "\t E Emission: " << CV::E_emis) << "\n\n";

	//double** XestS = new double* [nEnode];
	//double** XdecS = new double* [nEnode];
	//double** Xs = new double* [nEnode];
	EV::val_num_est = new double[nPlt]();
	EV::val_num_decom = new double[nPlt]();
	EV::val_Xop = new double* [nEnode];
	EV::val_Xest = new double* [nEnode];
	EV::val_Xdec = new double* [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_Xop[n] = new double[nPlt]();
		EV::val_Xest[n] = new double[nPlt]();
		EV::val_Xdec[n] = new double[nPlt]();
		//XestS[n] = new double[nPlt]();
		//XdecS[n] = new double[nPlt]();
		//Xs[n] = new double[nPlt]();
		for (int i = 0; i < nPlt; i++)
		{
			EV::val_num_est[i] += EV::Xest[n][i].get(GRB_DoubleAttr_X);
			EV::val_num_decom[i] += EV::Xdec[n][i].get(GRB_DoubleAttr_X);
			EV::val_Xop[n][i] = std::round(EV::Xop[n][i].get(GRB_DoubleAttr_X));
			EV::val_Xest[n][i] = std::round(EV::Xest[n][i].get(GRB_DoubleAttr_X));
			EV::val_Xdec[n][i] = std::round(EV::Xdec[n][i].get(GRB_DoubleAttr_X));
			/*if (EV::val_Xest[n][i] > 0)
			{
				cout << "Xest[" << n << "][" << i << "] = " << EV::val_Xest[n][i] << endl;
			}
			if (EV::val_Xop[n][i] > 0)
			{
				cout << "Xop[" << n << "][" << i << "] = " << EV::val_Xop[n][i] << endl;
			}
			if (EV::val_Xdec[n][i] > 0)
			{
				cout << "Xdec[" << n << "][" << i << "] = " << EV::val_Xdec[n][i] << endl;
			}*/

			//XestS[n][i] = EV::Xest[n][i]);
			//XdecS[n][i] = EV::Xdec[n][i]);
			//Xs[n][i] = EV::Xop[n][i]);
			//if (XestS[n][i] > 0)
			//{
			//	//std::cout << "Xest[" << n << "][" << i << "] = " << XestS[n][i] << endl;
			//	fid << "Xest[" << n << "][" << i << "] = " << XestS[n][i] << endl;
			//}
			//if (XdecS[n][i] > 0)
			//{
			//	//std::cout << "Xdec[" << n << "][" << i << "] = " << XdecS[n][i] << endl;
			//	fid << "Xdec[" << n << "][" << i << "] = " << XdecS[n][i] << endl;
			//}
			//if (Xs[n][i] > 0)
			//{
			//	//cout << "X[" << n << "][" << i << "] = " << Xs[n][i] << endl;
			//	fid << "Xop[" << n << "][" << i << "] = " << Xs[n][i] << endl;
			//}
		}
	}

	EV::val_flowE = new double* [nBr];
	for (int b = 0; b < nBr; b++)
	{
		EV::val_flowE[b] = new double[Te.size()]();
		for (int t = 0; t < Te.size(); t++)
		{
			//if (t > periods2print) { break; }
			EV::val_total_flow += time_weight[t] * EV::flowE[b][t].get(GRB_DoubleAttr_X);
			EV::val_flowE[b][t] = EV::flowE[b][t].get(GRB_DoubleAttr_X);
			//if (std::abs(fs[b][t]) > 10e-3)
			//{
			//	//std::cout << "flowE[" << n << "][" << t << "][" << Enodes[n].adj_buses[m] << "] = " << fs[n][t][m] << endl;
			//	cout << "flowE[" << b << "][" << t << "] = " << fs[b][t] << endl;
			//}
		}
	}

	EV::val_Ze = new double[nBr]();
	EV::val_num_est_trans = 0;
	for (int b = 0; b < nBr; b++)
	{
		EV::val_num_est_trans += EV::Ze[b].get(GRB_DoubleAttr_X);
		EV::val_Ze[b] = std::round(EV::Ze[b].get(GRB_DoubleAttr_X));

		//if (EV::val_Ze[b] > 10e-3)
		//{
		//	//std::cout << "Z[" << n << "][" << Enodes[n].adj_buses[m] << "] = " << Zs[n][m] << endl;
		//	std::cout << "Z[" << b << "] = " << EV::val_Ze[b] << endl;
		//}
	}

	EV::val_prod = new double** [nEnode];
	EV::val_X = new double** [nEnode];
	double total_yearly_prod = 0;
	EV::val_total_prod = new double[nPlt]();
	double flowGE = 0;
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_prod[n] = new double* [Te.size()];
		EV::val_X[n] = new double* [Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			//if (t > periods2print) { break; }

			EV::val_prod[n][t] = new double[nPlt]();
			EV::val_X[n][t] = new double[nPlt]();
			for (int i = 0; i < nPlt; i++)
			{
				/*
				if (Plants[i].type == "ng" || Plants[i].type == "CT" || Plants[i].type == "CC" || Plants[i].type == "CC-CCS")
				{
					double s1 = time_weight[t] * EV::prod[n][t][i].get(GRB_DoubleAttr_X);
					flowGE += s1 * Plants[i].heat_rate;
				}*/
				EV::val_total_prod[i] += time_weight[t] * EV::prod[n][t][i].get(GRB_DoubleAttr_X);
				total_yearly_prod += time_weight[t] * EV::prod[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_prod[n][t][i] = EV::prod[n][t][i].get(GRB_DoubleAttr_X);
				EV::val_X[n][t][i] = EV::X[n][t][i].get(GRB_DoubleAttr_X);
				//if (prodS[n][t][i] > 10e-3)
				//{
				//	//	std::cout << "prod[" << n << "][" << t << "][" << i << "] = " << prodS[n][t][i] << endl;
				//	fid << "prod[" << n << "][" << t << "][" << i << "] = " << prodS[n][t][i] << endl;
				//}
			}
		}
	}
	for (int i = 0; i < nPlt; i++)
	{
		if (i == 14)
		{
			int gg = 0;
		}
		if (EV::val_total_prod[i] > 10)
		{
			cout << Plants[i].type << ": \t " << EV::val_total_prod[i] << endl;
		}
	}


	cout << "\t\t total yearly product: " << total_yearly_prod << endl;
	EV::val_curtE = new double* [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		//cout << EV::Xop[n][12].get(GRB_DoubleAttr_X)<<endl; // wind offshore new: num of operating units
		EV::val_curtE[n] = new double[Te.size()]();
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_total_curt += time_weight[t] * EV::curtE[n][t].get(GRB_DoubleAttr_X);

			EV::val_curtE[n][t] = EV::curtE[n][t].get(GRB_DoubleAttr_X);
			//if (curtS[n][t] > 10e-3)
			//{
			//	//	std::cout << "curtE[" << n << "][" << t << "] = " << curtS[n][t] << endl;
			//	fid << "curtE[" << n << "][" << t << "] = " << curtS[n][t] << endl;
			//}
		}
	}

	EV::val_theta = new double* [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_theta[n] = new double[Te.size()]();
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_theta[n][t] = EV::theta[n][t].get(GRB_DoubleAttr_X);
		}
	}


	// Storage variables
	double** YeCDs = new double* [nEnode];
	EV::val_storage_cap = 0;
	for (int n = 0; n < nEnode; n++)
	{
		YeCDs[n] = new double[neSt]();
		for (int r = 0; r < neSt; r++)
		{
			EV::val_storage_cap += EV::YeCD[n][r].get(GRB_DoubleAttr_X);
			YeCDs[n][r] = EV::YeCD[n][r].get(GRB_DoubleAttr_X);
			/*if (YeCDs[n][r] > 10e-3)
			{
				fid << "Ye_ch[" << n << "][" << r << "] = " << YeCDs[n][r] << endl;
			}*/
		}
	}
	double** YeLevs = new double* [nEnode];
	EV::val_storage_lev = 0;
	for (int n = 0; n < nEnode; n++)
	{
		YeLevs[n] = new double[neSt]();
		for (int r = 0; r < neSt; r++)
		{
			EV::val_storage_lev += EV::YeLev[n][r].get(GRB_DoubleAttr_X);
			YeLevs[n][r] = EV::YeLev[n][r].get(GRB_DoubleAttr_X);
			/*if (YeLevs[n][r] > 10e-3)
			{
				std::cout << "Ye_lev[" << n << "][" << r << "] = " << YeLevs[n][r] << endl;
			}*/
		}
	}

	EV::val_YeStr = new double* [nEnode];
	EV::val_num_storage = 0;
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_YeStr[n] = new double[neSt]();
		for (int r = 0; r < neSt; r++)
		{
			EV::val_num_storage += EV::YeStr[n][r].get(GRB_DoubleAttr_X);
			EV::val_YeStr[n][r] = EV::YeStr[n][r].get(GRB_DoubleAttr_X);
			//YeStr[n][r] = EV::YeStr[n][r]);
			/*if (EV::val_YeStr[n][r] > 10e-3)
			{
				std::cout << "Ye_Str[" << n << "][" << r << "] = " << EV::val_YeStr[n][r] << endl;
			}*/
		}
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
				//if (eSchS[n][t][i] > 10e-3)
				//{
				//	//	std::cout << "prod[" << n << "][" << t << "][" << i << "] = " << prodS[n][t][i] << endl;
				//	fid << "eS_ch[" << n << "][" << t << "][" << i << "] = " << eSchS[n][t][i] << endl;
				//}
			}
		}
	}

	EV::val_eSlev = new double** [nEnode];
	for (int n = 0; n < nEnode; n++)
	{
		EV::val_eSlev[n] = new double* [Te.size()];
		for (int t = 0; t < Te.size(); t++)
		{
			EV::val_eSlev[n][t] = new double[neSt]();
			for (int i = 0; i < neSt; i++)
			{
				EV::val_eSlev[n][t][i] = EV::eSlev[n][t][i].get(GRB_DoubleAttr_X);
			}
		}
	}

	//for (int n = 0; n < nEnode; n++)
	//{
	//	for (int t = 0; t < Te.size(); t++)
	//	{
	//		double s1 = EV::theta[n][t].get(GRB_DoubleAttr_X);
	//		if (s1>0)
	//		{
	//			//cout << "theta[" << n << "][" << t << "] = " << s1 << endl;
	//		}
	//	}
	//}
	//fid.close();
#pragma endregion
}


void Get_GV_vals(GRBModel Model)
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
	double pi = 3.1415;
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

#pragma region print NG network variables
	/*ofstream fid2;
	if (Setting::DGSP_active)
	{
		fid2.open("DGSP.txt");
	}
	else
	{
		fid2.open("DVg.txt");
	}*/

	GV::val_NG_system_cost = GV::NG_system_cost.get(GRB_DoubleAttr_X);
	GV::val_pipe_cost = GV::pipe_cost.get(GRB_DoubleAttr_X);
	GV::val_NG_import_cost = GV::NG_import_cost.get(GRB_DoubleAttr_X);
	GV::val_strInv_cost = GV::strInv_cost.get(GRB_DoubleAttr_X);
	GV::val_gStrFOM_cost = GV::gStrFOM_cost.get(GRB_DoubleAttr_X);
	GV::val_ngShedd_cost = GV::gShedd_cost.get(GRB_DoubleAttr_X);
	GV::val_rngShedd_cost = GV::rngShedd_cost.get(GRB_DoubleAttr_X);

	CV::val_NG_emis = CV::NG_emis.get(GRB_DoubleAttr_X);


	/*fid2 << "Elapsed time: " << Elapsed_time << endl;
	fid2 << "\t Total cost for both networks:" << obj << endl;
	fid2 << "\t NG Network Obj Value:" << GV::NG_system_cost) << endl;
	fid2 << "\t Gap: " << gap << " Status:" << cplex.getStatus() << endl;
	fid2 << "\t NG import Cost: " << GV::NG_import_cost) << endl;
	fid2 << "\t Pipeline Establishment Cost: " << GV::pipe_cost) << endl;
	fid2 << "\t Storatge Investment Cost: " << GV::strInv_cost) << endl;
	fid2 << "\t NG Storage Cost: " << GV::gStrFOM_cost) << endl;
	fid2 << "\t NG Load Shedding Cost: " << GV::gShedd_cost) << endl;
	fid2 << "\t NG Emission: " << CV::NG_emis) << endl;*/
	//fid2 << endl;
	GV::val_supply = new double* [nGnode];
	GV::val_total_nodal_supply = new double[nGnode]();
	for (int k = 0; k < nGnode; k++)
	{
		GV::val_supply[k] = new double[Tg.size()]();
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			GV::val_total_nodal_supply[k] += RepDaysCount[tau] * GV::supply[k][tau].get(GRB_DoubleAttr_X);
			GV::val_supply[k][tau] = GV::supply[k][tau].get(GRB_DoubleAttr_X);
			/*if (s1 > 0)
			{
				fid2 << "supply[" << k << "][" << tau << "] = " << s1 << endl;
			}*/
		}
	}
	//	fid2 << endl;
	//GV::val_g_curt = new double[nGnode]();
	GV::val_ng_curt = 0; GV::val_rng_curt = 0;
	for (int k = 0; k < nGnode; k++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			GV::val_ng_curt += RepDaysCount[tau] * GV::curtG[k][tau].get(GRB_DoubleAttr_X);
			GV::val_rng_curt += RepDaysCount[tau] * GV::curtRNG[k][tau].get(GRB_DoubleAttr_X);
			/*if (s1 > 0)
			{
				fid2 << "CurtG[" << k << "][" << tau << "] = " << s1 << endl;
			}*/
		}
	}
	//fid2 << endl;
	GV::val_Zg = new double[nPipe](); GV::val_num_est_pipe = 0;
	for (int i = 0; i < nPipe; i++)
	{
		GV::val_num_est_pipe += GV::Zg[i].get(GRB_DoubleAttr_X);
		GV::val_Zg[i] = GV::Zg[i].get(GRB_DoubleAttr_X);
		/*if (s1 > 0.01)
		{
			fid2 << "Zg[" << PipeLines[i].from_node << "][" << PipeLines[i].to_node << "] = " << s1 << endl;
		}*/
	}
	//fid2 << endl;
	GV::val_total_flowGG = 0;
	for (int i = 0; i < nPipe; i++)
	{
		for (int tau = 0; tau < Tg.size(); tau++)
		{
			GV::val_total_flowGG += RepDaysCount[tau] * GV::flowGG[i][tau].get(GRB_DoubleAttr_X);
			/*if (s1 > 0.001)
			{
				fid2 << "flowGG[" << PipeLines[i].from_node << "][" << PipeLines[i].to_node << "][" << tau << "]= " << s1 << endl;
			}*/
		}
	}
	//fid2 << endl;
	GV::val_flowGE = new double** [nGnode];
	GV::val_total_flowGE = 0;
	for (int k = 0; k < nGnode; k++)
	{
		GV::val_flowGE[k] = new double* [nEnode];
		for (int kp : Gnodes[k].adjE)
		{
			GV::val_flowGE[k][kp] = new double[Tg.size()]();
			for (int tau = 0; tau < Tg.size(); tau++)
			{
				GV::val_flowGE[k][kp][tau] = GV::flowGE[k][kp][tau].get(GRB_DoubleAttr_X);
				GV::val_total_flowGE += RepDaysCount[tau] * GV::flowGE[k][kp][tau].get(GRB_DoubleAttr_X);
				/*if (std::abs(GV::val_flowGE[k][kp][tau]) > 0.1)
				{
					std::cout << "flowGE[" << k << "][" << kp << "][" << tau << "] = " << GV::val_flowGE[k][kp][tau] << endl;
				}*/
			}
		}
	}

	////fid2 << endl;
	GV::val_total_flowGL = 0;
	for (int k = 0; k < nGnode; k++)
	{
		for (int kp : Gnodes[k].adjS)
		{
			for (int tau = 0; tau < Tg.size(); tau++)
			{
				GV::val_total_flowGL += RepDaysCount[tau] * GV::flowGL[k][kp][tau].get(GRB_DoubleAttr_X);
				/*if (std::abs(s1) > 0.1)
				{
					fid2 << "flowGL[" << k << "][" << kp << "][" << tau << "] = " << s1 << endl;
				}*/
			}
		}
	}
	//fid2 << endl;
	GV::val_total_flowVG = 0;
	for (int k = 0; k < nGnode; k++)
	{
		for (int kp : Gnodes[k].adjS)
		{
			for (int tau = 0; tau < Tg.size(); tau++)
			{
				GV::val_total_flowVG += RepDaysCount[tau] * GV::flowVG[kp][k][tau].get(GRB_DoubleAttr_X);
				/*if (std::abs(s1) > 0.1)
				{
					fid2 << "flowVG[" << kp << "][" << k << "][" << tau << "] = " << s1 << endl;
				}*/
			}
		}
	}

	//fid2 << endl;
	GV::val_storage = new double[nSVL]();
	GV::val_Xstr = new double[nSVL]();
	for (int j = 0; j < nSVL; j++)
	{
		GV::val_storage[j] += GV::Xstr[j].get(GRB_DoubleAttr_X);
		GV::val_Zg[j] = GV::Xstr[j].get(GRB_DoubleAttr_X);
		/*if (s1 > 0.001)
		{
			fid2 << "Xstr[" << j << "]= " << s1 << endl;
		}*/
	}
	int gg = 0;
	//fid2 << endl;
	GV::val_vapor = new double[nSVL]();
	for (int j = 0; j < nSVL; j++)
	{
		GV::val_vapor[j] += GV::Xvpr[j].get(GRB_DoubleAttr_X);

		/*if (s1 > 0.001)
		{
			fid2 << "Xvpr[" << j << "]= " << s1 << endl;
		}*/
	}


	//fid2 << endl;
	//for (int j = 0; j < nSVL; j++)
	//{
	//	for (int t = 0; t < Tg.size(); t++)
	//	{
	//		double s1 = GV::Sstr[j][t]);
	//		if (s1 > 0.001)
	//		{
	//			fid2 << "Sstr[" << j << "][" << t << "]= " << s1 << endl;
	//		}
	//	}
	//}
	//fid2 << endl;
	//for (int j = 0; j < nSVL; j++)
	//{
	//	for (int t = 0; t < Tg.size(); t++)
	//	{
	//		double s1 = GV::Svpr[j][t]);
	//		if (s1 > 0.001)
	//		{
	//			fid2 << "Svpr[" << j << "][" << t << "]= " << s1 << endl;
	//		}
	//	}
	//}
	//fid2 << endl;
	//for (int j = 0; j < nSVL; j++)
	//{
	//	for (int t = 0; t < Tg.size(); t++)
	//	{
	//		double s1 = GV::Sliq[j][t]);
	//		if (s1 > 0.001)
	//		{
	//			fid2 << "Sliq[" << j << "][" << t << "]= " << s1 << endl;
	//		}
	//	}
	//}
	//fid2.close();
#pragma endregion

}



void Print_Results(double Elapsed_time, double status)
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

#pragma region Print all decision varialbes of both networks
	ofstream fid0;
	// only for approach 1 or  case 1
	int is_given = (std::round(Setting::is_xi_given));
	int em_lim = 100 * (Setting::Emis_redu_goal);
	int rps = 100 * (Setting::RPS);
	int rng = 100 * Setting::RNG_cap;
	/*string name = "Rep " + std::to_string(Setting::Num_rep_days) +
		"-A1" + "-case " + std::to_string(Setting::Case) + "-xi_given " + std::to_string(is_given)
		+ "-emis_lim " + std::to_string(em_lim) + "-RPS " + std::to_string(rps)
		+ "-RNG " + std::to_string(rng) + ".csv";*/
	string name = "Rep " + std::to_string(Setting::Num_rep_days) +
		+ "-case " + std::to_string(Setting::Case) 
		+ "-emis_lim " + std::to_string(em_lim) + "-RPS " + std::to_string(rps)
		+ "-RNG " + std::to_string(rng) + ".csv";

	if (Setting::print_all_vars)
	{
		fid0.open(name, std::ios::app);

		for (int i = 0; i < nPlt; i++)
		{
			fid0 << "\nProd[" << i << "] " << ",";

			for (int t = 0; t < Te.size(); t++)
			{
				double node_prod = 0;
				for (int n = 0; n < nEnode; n++)
				{
					node_prod += EV::val_prod[n][t][i];
				}
				fid0 << node_prod << ",";
			}
		}


		fid0 << "\n\nBattery Charge" << ",";

		for (int t = 0; t < Te.size(); t++)
		{
			double all_node_charge = 0;
			for (int n = 0; n < nEnode; n++)
			{
				all_node_charge += EV::val_eSch[n][t][0];
			}
			fid0 << all_node_charge << ",";
		}


		fid0 << "\nBattery Discharge" << ",";

		for (int t = 0; t < Te.size(); t++)
		{
			double all_node_discharge = 0;
			for (int n = 0; n < nEnode; n++)
			{
				all_node_discharge += EV::val_eSdis[n][t][0];
			}
			fid0 << all_node_discharge << ",";
		}


		fid0 << "\n\nCurt";
		for (int t = 0; t < Te.size(); t++)
		{
			double curT = 0;
			for (int n = 0; n < nEnode; n++)
			{
				curT += EV::val_curtE[n][t];
			}
			fid0 << "," << curT;
		}

		for (int i = 0; i < neSt; i++)
		{
			fid0 << "\n eSlev[" << i << "] " << ",";

			for (int t = 0; t < Te.size(); t++)
			{
				double node_prod = 0;
				for (int n = 0; n < nEnode; n++)
				{
					node_prod += EV::val_eSlev[n][t][i];
				}
				fid0 << node_prod << ",";
			}
		}
		fid0.close();
	}


#pragma endregion


#pragma region Print in CSV file
	ofstream fid;
	name = "NGES_Results.csv";
	fid.open(name, std::ios::app);
	if (Setting::print_results_header)
	{
		// problem setting
		fid << "\nSol_Time" << ",";
		fid << "Num_Rep_days" << ",";
		fid << "Approach :" << ",";
		fid << "Case: " << ",";
		fid << "xi_given" << ",";
		fid << "xi_val" << ",";
		fid << "Emis_lim" << ",";
		fid << "RPS" << ",";
		fid << "RNG_cap" << ",";
		fid << "MIO Gap:" << ",";
		fid << "Benders" << ",";
		fid << "Benders Multicut" << ",";
		fid << "Benders Iter" << ",";
		fid << "UC Active" << ",";
		fid << "UC Vars Relaxed" << ",";


		// electricity
		fid << "Total_E_cost" << ",";
		fid << "Est_plant_cost" << ",";
		fid << "Est_trans_cost" << ",";
		fid << "Decom_cost" << ",";
		fid << "Fixed_cost" << ",";
		fid << "Startup_cost" << ",";
		fid << "Variable_cost" << ",";
		fid << "dfo_coal_or_NG_emis_cost" << ",";
		fid << "Coal_dfo_or_NG_fuel_cost" << ",";
		fid << "Shedding_cost" << ",";
		fid << "Storage_cost" << ",";
		fid << "Num_storage" << ",";
		fid << "Storage_level" << ",";
		fid << "Storage_cap" << ",";
		fid << "Total_curt" << ",";
		fid << "Num_new_trans" << ",";
		fid << "Total_flow" << ",";
		fid << "E_emis" << ",,";

		// natural gas
		fid << "Total_NG_cost" << ",";
		fid << "Import_cost" << ",";
		fid << "Pipe_est_cost" << ",";
		fid << "Str_est_cost" << ",";
		fid << "Str_fix_cost" << ",";
		fid << "NG_shedding_cost" << ",";
		fid << "RNG_shedding_cost" << ",";
		fid << "Num_est_pipe" << ",";
		fid << "Total_NG_shedding" << ",";
		fid << "Total_RNG_shedding" << ",";
		fid << "flowGG" << ",";
		fid << "flowGE" << ",";
		fid << "flowGL" << ",";
		fid << "flowVG" << ",";
		fid << "NG_emis" << ",";
		// vector variables of both networks
		fid << ",,";
		fid << "num_est(nplt)" << ","; // vector of nplt
		fid << "num_decom(nplt)" << ","; // vector of nplt
		fid << "prod(nplt)" << ","; // vector of nplt

		fid << "supply(nGnode)" << ",";
		fid << "storage(nSLV)" << ",";
		fid << "vapor(nSLV)" << ",";


	}
	fid << "\n" << Elapsed_time << ",";
	fid << Setting::Num_rep_days << ",";
	if (Setting::Approach_1_active)
	{
		fid << 1 << ",";
	}
	else if (Setting::Approach_2_active)
	{
		fid << 2 << ",";
	}
	else
	{
		fid << "decoup" << ",";
	}
	fid << Setting::Case << ",";
	//fid << Setting::Approach_1_active << ",";
	//fid << Setting::Approach_2_active << ",";

	fid << Setting::is_xi_given << ",";
	fid << CV::val_xi << ",";
	fid << Setting::Emis_redu_goal << ",";
	fid << Setting::RPS << ",";
	fid << Setting::RNG_cap << ",";
	//if (status == -1)// problem is infeasible
	//{
	//	fid.close();
	//	return;
	//}
	fid << EV::MIP_gap << ",";
	fid << Setting::use_benders << ",";
	fid << Setting::multi_cut_active << ",";
	fid << EV::Benders_iter << ",";
	fid << Setting::UC_active << ",";
	fid << Setting::relax_UC_vars << ",";
	fid << EV::val_e_system_cost << ",";
	fid << EV::val_est_cost << ",";
	fid << EV::val_est_trans_cost << ",";
	fid << EV::val_decom_cost << ",";
	fid << EV::val_fixed_cost << ",";
	fid << EV::val_startup_cost << ",";
	fid << EV::val_var_cost << ",";
	fid << EV::val_dfo_coal_emis_cost << ",";
	fid << EV::val_thermal_fuel_cost << ",";
	fid << EV::val_shedding_cost << ",";
	fid << EV::val_elec_storage_cost << ",";
	fid << EV::val_num_storage << ",";
	fid << EV::val_storage_lev << ",";
	fid << EV::val_storage_cap << ",";
	fid << EV::val_total_curt << ",";
	fid << EV::val_num_est_trans << ",";
	fid << EV::val_total_flow << ",";
	fid << CV::val_E_emis << ",,";

	// NG
	fid << GV::val_NG_system_cost << ",";
	fid << GV::val_NG_import_cost << ",";
	fid << GV::val_pipe_cost << ",";
	fid << GV::val_strInv_cost << ",";
	fid << GV::val_gStrFOM_cost << ",";
	fid << GV::val_ngShedd_cost << ",";
	fid << GV::val_rngShedd_cost << ",";
	fid << GV::val_num_est_pipe << ",";
	fid << GV::val_ng_curt << ",";
	fid << GV::val_rng_curt << ",";
	fid << GV::val_total_flowGG << ",";
	fid << GV::val_total_flowGE << ",";
	fid << GV::val_total_flowGL << ",";
	fid << GV::val_total_flowVG << ",";
	fid << CV::val_NG_emis << ",";

	fid << ",";
	for (int i = 0; i < Params::Plants.size(); i++)
	{
		fid << "," << EV::val_num_est[i];
	}
	fid << ",";
	for (int i = 0; i < Params::Plants.size(); i++)
	{
		fid << "," << EV::val_num_decom[i];
	}
	fid << ",";
	for (int i = 0; i < Params::Plants.size(); i++)
	{
		fid << "," << EV::val_total_prod[i];
	}
	fid << ",";
	for (int j = 0; j < Params::Gnodes.size(); j++)
	{
		fid << "," << GV::val_total_nodal_supply[j];
	}
	fid << ",";
	for (int j = 0; j < Params::SVLs.size(); j++)
	{
		fid << "," << GV::val_storage[j];
	}
	fid << ",";
	for (int j = 0; j < Params::SVLs.size(); j++)
	{
		fid << "," << GV::val_vapor[j];
	}
	fid.close();
#pragma endregion

#pragma region Print in a text file
	ofstream fid2;
	fid2.open("NGES_Results.txt", std::ios::app);
	// problem setting
	fid2 << "\n ************************************************************************* ";
	fid2 << "\n Num_Rep_days: " << Setting::Num_rep_days;
	fid2 << "\t Approach 1: " << Setting::Approach_1_active;
	fid2 << "\tApproach 2: " << Setting::Approach_2_active;
	fid2 << "\t Case: " << Setting::Case;
	fid2 << "\txi_given: " << Setting::is_xi_given;
	fid2 << "\txi_val: " << CV::val_xi;
	fid2 << "\tEmis_lim: " << Setting::Emis_redu_goal;
	fid2 << "\tRPS: " << Setting::RPS;
	fid2 << "\tRNG_cap: " << Setting::RNG_cap;
	fid2 << "\tBenders: " << Setting::use_benders;
	fid2 << "\tBenders Multicut: " << Setting::multi_cut_active;
	fid2 << "\tUC Active: " << Setting::UC_active;
	fid2 << "\tUC vars relaxed: " << Setting::relax_UC_vars;


	fid2 << "\n\t Elapsed time: " << Elapsed_time << endl;
	fid2 << "\t Total cost for both networks:" << EV::val_e_system_cost + GV::val_NG_system_cost << endl;

	fid2 << "\n \t Electricity Network Obj Value: " << EV::val_e_system_cost;
	fid2 << "\n \t Plant Establishment Cost: " << EV::val_est_cost;
	fid2 << "\n \t Trnas. Establishment Cost: " << EV::val_est_trans_cost;
	fid2 << "\n \t Decommissioning Cost: " << EV::val_decom_cost;
	fid2 << "\n \t Fixed Cost: " << EV::val_fixed_cost;
	fid2 << "\n \t Startup Cost: " << EV::val_startup_cost;
	fid2 << "\n \t Variable Cost: " << EV::val_var_cost;
	fid2 << "\n \t [NG] emission cost: " << EV::val_dfo_coal_emis_cost;
	fid2 << "\n \t ([NG and] nuclear) Fuel Cost: " << EV::val_thermal_fuel_cost;
	fid2 << "\n \t Load Shedding Cost: " << EV::val_shedding_cost;
	fid2 << "\n \t Storage Cost: " << EV::val_elec_storage_cost;
	fid2 << "\n \t E Emission: " << CV::val_E_emis << "\n";


	fid2 << "\n\n\t NG Network Obj Value:" << GV::val_NG_system_cost << endl;
	fid2 << "\t NG import Cost: " << GV::val_NG_import_cost << endl;
	fid2 << "\t Pipeline Establishment Cost: " << GV::val_pipe_cost << endl;
	fid2 << "\t Storatge Investment Cost: " << GV::val_strInv_cost << endl;
	fid2 << "\t NG Storage Cost: " << GV::val_gStrFOM_cost << endl;
	fid2 << "\t NG Load NG Shedding Cost: " << GV::val_ngShedd_cost << endl;
	fid2 << "\t NG Load RNG Shedding Cost: " << GV::val_rngShedd_cost << endl;
	fid2 << "\t NG Emission: " << CV::val_NG_emis << endl;
	fid2 << "\n\n\n";
	fid2.close();
#pragma endregion

}


