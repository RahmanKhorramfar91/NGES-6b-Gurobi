#include "ProblemData.h"
vector<gnode> Params::Gnodes;
vector<pipe> Params::PipeLines;
vector<enode> Params::Enodes;
vector<eStore> Params::Estorage;
vector<exist_gSVL> Params::Exist_SVL;
vector<SVL> Params::SVLs;
vector<plant> Params::Plants;
vector<branch> Params::Branches;
vector<int> Params::Tg;
vector<int> Params::Te;
vector<int> Params::RepDaysCount;
vector<int> Params::time_weight;
double Params::WACC;
int Params::trans_unit_cost;
int Params::trans_line_lifespan;
double Params::NG_price;
double Params::dfo_pric;
double Params::coal_price;
double Params::nuclear_price;
double Params::E_curt_cost;
double Params::G_curt_cost;
double Params::RNG_cap;
double Params::RNG_price;
double Params::pipe_per_mile;
double Params::NG_emis_rate;
int Params::pipe_lifespan;
int Params::SVL_lifetime;
int Params::battery_lifetime;
int Params::Num_Rep_Days;
map<int, vector<int>> Params::Le;
map<int, vector<int>> Params::Lg;


void Read_rep_days(string name, vector<int>& Rep, vector<int>& RepCount)
{
	ifstream fid(name);
	string line;
	while (getline(fid, line))
	{
		double re, rc;
		std::istringstream iss(line);
		iss >> re >> rc;
		Rep.push_back((int)re);
		RepCount.push_back((int)rc);
	}
}

vector<enode> enode::read_bus_data(string name)
{
	vector<enode> Buses;


	ifstream fid(name);
	string line;
	while (getline(fid, line))
	{
		std::istringstream iss(line);
		double bus_num;
		iss >> bus_num;
		enode bs((int)bus_num);
		Buses.push_back(bs);
	}

	fid.close();
	return Buses;
}

void  enode::read_adj_data(string name, vector<enode>& Enodes)
{
	ifstream fid(name);
	string line;
	int j = -1;
	while (getline(fid, line))
	{
		string temp = "";
		j++;
		int if_two_tab = 0;
		for (int i = 0; i < line.length(); i++)
		{
			if (line[i] == '\t' || i == line.length() - 1)
			{
				if_two_tab++;
				if (if_two_tab > 1)
				{
					break;
				}
				Enodes[j].adj_buses.push_back(std::stod(temp));

				temp = "";
			}
			else
			{
				temp.push_back(line[i]);
				if_two_tab = 0;
			}
		}
	}


}

void enode::read_exist_plt_data(string FileName, vector<enode>& Enodes)
{

	std::map<string, int> sym2ind = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7} };


	ifstream fid(FileName);
	string line;
	while (getline(fid, line))
	{
		std::istringstream iss(line);
		string type;
		double bus_num, Pmax, Pmin, g1, g2, g3, count;
		iss >> bus_num >> type >> Pmax >> Pmin >> g1 >> g2 >> g3 >> count;
		int bn = (int)bus_num;
		Enodes[bn].Init_plt_types.push_back(type);
		int ind1 = sym2ind[type];
		Enodes[bn].Init_plt_ind.push_back(ind1);
		double* plim = new double[2]{ Pmin,Pmax };
		Enodes[bn].Init_plt_prod_lim.push_back(plim);
		Enodes[bn].Init_plt_count.push_back((int)count);
	}
}

void enode::read_demand_data(string name, vector<enode>& Enodes)
{
	ifstream fid(name);
	string line;
	while (getline(fid, line))
	{
		//vector<string> v;
		string temp = "";
		int j = 0;
		for (int i = 0; i < line.length(); i++)
		{
			if (line[i] == '\t' || i == line.length() - 1)
			{
				Enodes[j].demand.push_back(std::stof(temp));
				j++;
				temp = "";
			}
			else
			{
				temp.push_back(line[i]);
			}
		}
	}
	fid.close();
}


vector<plant> plant::read_new_plant_data(string name)
{
	vector<plant> NewPlants;
	ifstream fid(name);
	string line;
	while (getline(fid, line))
	{
		//	string t, int n, int ise, int cap, int f, int v, double emi, double hr, int lt, int dec,
			//	double pmax, double ru, double rd, int emic
		std::istringstream iss(line);
		double n, ise, cap, f, v, emi, hr, lt, dec, pmax, pmin, ru, rd, emic, numM;
		//double  n, capex,pmax, pmin,ru,rd, fix_cost, var_cost, decom_cost, emis_cost, lifespan;
		string type;
		double h, emis_rate;
		iss >> type >> n >> ise >> cap >> f >> v >> emi >> hr >> lt >> dec >> pmax >> pmin >> ru >> rd >> emic >> numM;
		//iss >> type >> n >> capex >>pmax>>pmin>>ru>>rd>> fix_cost >> var_cost >> h >> emis_rate >> decom_cost >> emis_cost >> lifespan;
		plant np(type, (int)n, (int)ise, (int)cap, (int)f, (int)v, emi, hr, (int)lt, (int)dec, pmax, pmin, ru, rd, (int)emic, (int)numM);
		NewPlants.push_back(np);
	}
	fid.close();
	return NewPlants;
}

vector<eStore> eStore::read_elec_storage_data(string name)
{
	vector<eStore> Estorage;
	ifstream fid(name);
	string line;
	while (getline(fid, line))
	{
		std::istringstream iss(line);
		double en, pow, ch, dis, fom;
		iss >> en >> pow >> ch >> dis >> fom;
		eStore str((int)en, (int)pow, ch, dis, fom);
		Estorage.push_back(str);
	}
	fid.close();
	return Estorage;
}

void plant::read_VRE_profile(string FN1, string FN2, string FN3, vector<plant>& Plants)
{
	std::map<string, int> sym2pltType = { {"ng",0},{"dfo", 1},
{"solar", 2},{"wind", 3},{"wind_offshore", 4},{"hydro", 5},{"coal",6},{"nuclear",7},
		{"Ct",8},{"CC",9},{"CC-CCS",10},{"solar-UPV",11},{"wind-new",12},
		{"wind-offshore-new",13},{"hydro-new",14},{"nuclear-new",15} };

	ifstream fid1("profile_hydro_hourly.txt");
	ifstream fid2("profile_wind_hourly.txt");
	ifstream fid3("profile_solar_hourly.txt");
	string line;
	while (getline(fid1, line))
	{
		std::istringstream iss(line);
		double pp;
		iss >> pp;
		Plants[sym2pltType["hydro"]].prod_profile.push_back(pp);
		Plants[sym2pltType["hydro-new"]].prod_profile.push_back(pp);
	}

	while (getline(fid2, line))
	{
		std::istringstream iss(line);
		double pp;
		iss >> pp;
		Plants[sym2pltType["wind"]].prod_profile.push_back(pp);
		Plants[sym2pltType["wind-new"]].prod_profile.push_back(pp);
	}
	while (getline(fid3, line))
	{
		std::istringstream iss(line);
		double pp;
		iss >> pp;
		Plants[sym2pltType["solar"]].prod_profile.push_back(pp);
		Plants[sym2pltType["solar-UPV"]].prod_profile.push_back(pp);
	}

}

vector<branch> branch::read_branch_data(int nBus, string FN, map<int, vector<int>>& Le)
{

	vector<branch> Branches;
	ifstream fid(FN); // to from, max flow, 

	string line;
	int lc = 0;
	while (getline(fid, line))
	{
		// to , from, is_exist, maxflow, susceptance, distance
		std::istringstream iss(line);
		double to, from, ise, maxF, sus, dist;
		iss >> from >> to >> ise >> maxF >> sus >> dist;
		branch nb(from, to, dist, ise, maxF, sus);
		//Lnm[i*200+ (int)s0] = vector<int>();
		Le[from * 200 + (int)to].push_back(lc);
		lc++;
		Branches.push_back(nb);
	}
	return Branches;

}


//vector<branch> branch::read_branch_data_old_version(int nBus, string FN,
//	string FN0, string FN1, string FN2, string FN3, string FN4, map<int, vector<int>>& Le)
//{
//
//	vector<branch> Branches;
//	ifstream fid(FN); // branch per node
//	ifstream fid0(FN0); // branch
//	ifstream fid1(FN1); // length
//	ifstream fid2(FN2); // is_exist
//	ifstream fid3(FN3); // max flow
//	ifstream fid4(FN4); // Susceptance
//
//	string line, line0, line1, line2, line3, line4;
//	int lc = 0;
//	for (int i = 0; i < nBus; i++)
//	{
//		getline(fid, line);// branch per node
//		getline(fid0, line0);// branch
//		getline(fid1, line1);// length
//		getline(fid2, line2); // is_exist
//		getline(fid3, line3); // max flow
//		getline(fid4, line4); // suscept
//
//		std::istringstream iss(line);
//		std::istringstream iss0(line0);
//		std::istringstream iss1(line1);
//		std::istringstream iss2(line2);
//		std::istringstream iss3(line3);
//		std::istringstream iss4(line4);
//
//		double sn;
//		iss >> sn;
//		if ((int)sn == 0) { continue; }
//
//
//		for (int j = 0; j < (int)sn; j++)
//		{
//			double s0, s1, s2, s3, s4;
//			iss0 >> s0;
//			iss1 >> s1;
//			iss2 >> s2;
//			iss3 >> s3;
//			iss4 >> s4;
//
//			branch nb(i, (int)s0, s1, s1, (int)s2, s3, s4);
//			//Lnm[i*200+ (int)s0] = vector<int>();
//			Le[i * 200 + (int)s0].push_back(lc);
//			lc++;
//			Branches.push_back(nb);
//		}
//	}
//	fid0.close();
//	fid1.close();
//	fid2.close();
//	fid3.close();
//
//	return Branches;
//
//}

vector<gnode> gnode::read_gnode_data(string name)
{
	vector<gnode>Gnodes;
	ifstream fid2(name);
	string line;
	while (getline(fid2, line))
	{
		double num, fips, outDem, capU, svl1, svl2;
		std::istringstream iss(line);

		iss >> num >> fips >> outDem >> capU >> svl1 >> svl2;
		vector<int> svls;
		svls.push_back((int)svl1);
		svls.push_back((int)svl2);
		gnode gn((int)num, (int)fips, (int)outDem, (int)capU, svls);
		Gnodes.push_back(gn);
	}
	fid2.close();
	return Gnodes;
}

vector<exist_gSVL> exist_gSVL::read_exist_SVL_data(string name)
{
	vector<exist_gSVL> Exist_SVL;
	ifstream fid;
	fid.open(name);
	string line;
	int count = 0;
	while (getline(fid, line))
	{
		double  a, b, c;
		std::istringstream iss(line);
		iss >> a >> b >> c;
		exist_gSVL esvl(count, (int)a, (int)b, (int)c);
		Exist_SVL.push_back(esvl);
		count++;
	}
	fid.close();
	return Exist_SVL;
}
vector<SVL> SVL::read_SVL_data(string name)
{

	vector<SVL> SVLs;
	ifstream fid;
	fid.open(name);
	string line;
	int count = 0;
	while (getline(fid, line))
	{
		double  a, b, c, d;
		std::istringstream iss(line);
		iss >> a >> b >> c >> d;
		SVL svl((int)a, (int)b, c, d);
		SVLs.push_back(svl);
		count++;
	}
	fid.close();
	return SVLs;
}
void gnode::read_Lexp_data(string FileName, vector<gnode>& Gnodes)
{
	ifstream fid2(FileName);
	string line;
	int j = -1;
	while (getline(fid2, line))
	{
		string temp = "";
		j++;
		int if_two_tab = 0;
		for (int i = 0; i < line.length(); i++)
		{
			if (line[i] == '\t' || i == line.length() - 1)
			{
				if_two_tab++;
				if (if_two_tab > 1) // end of number
				{
					break;
				}
				Gnodes[j].Lexp.push_back(std::stoi(temp));
				temp = "";
			}
			else
			{
				temp.push_back(line[i]);
				if_two_tab = 0;
			}
		}
	}

}
void gnode::read_Limp_data(string FileName, vector<gnode>& Gnodes)
{
	ifstream fid2(FileName);
	string line;
	int j = -1;
	while (getline(fid2, line))
	{
		if (line[0] == '\t')
		{
			j++;
			continue;
		}
		string temp = "";
		j++;
		int if_two_tab = 0;
		for (int i = 0; i < line.length(); i++)
		{
			if (line[i] == '\t' || i == line.length() - 1)
			{
				if_two_tab++;
				if (if_two_tab > 1) // end of number
				{
					break;
				}
				Gnodes[j].Limp.push_back(std::stoi(temp));
				temp = "";
			}
			else
			{
				temp.push_back(line[i]);
				if_two_tab = 0;
			}
		}
	}

}

void gnode::read_adjE_data(string FileName, vector<gnode>& Gnodes)
{
	ifstream fid2(FileName);
	string line;
	int j = -1;
	while (getline(fid2, line))
	{
		if (line == "\t")
		{
			j++;
			continue;
		}
		string temp = "";
		j++;
		int if_two_tab = 0;
		for (int i = 0; i < line.length(); i++)
		{
			if (line[i] == '\t' || i == line.length() - 1)
			{
				if_two_tab++;
				if (if_two_tab > 1) // end of number
				{
					break;
				}
				Gnodes[j].adjE.push_back(std::stoi(temp));
				temp = "";
			}
			else
			{
				temp.push_back(line[i]);
				if_two_tab = 0;
			}
		}
	}

}

vector<pipe> pipe::read_pipe_data(string name, map<int, vector<int>>& Lg)
{
	vector<pipe> PipeLines;
	ifstream fid2(name);
	string line;
	int lc = 0;
	while (getline(fid2, line))
	{
		double f, t, exi, le, ca;
		std::istringstream iss(line);
		iss >> f >> t >> exi >> le >> ca;
		pipe p2((int)f, (int)t, (int)exi, le, ca);
		Lg[(int)f * 200 + (int)t].push_back(lc);
		lc++;
		PipeLines.push_back(p2);
	}

	return PipeLines;
}


void gnode::read_ng_demand_data(string name, vector<gnode>& Gnodes)
{
	int nGnode = (int)Gnodes.size();

	ifstream fid2(name);
	string line;
	while (getline(fid2, line))
	{
		double dem2;
		std::istringstream iss(line);
		for (int i = 0; i < nGnode; i++)
		{
			iss >> dem2;
			Gnodes[i].demG.push_back(dem2);
		}
	}
}


void Read_Data()
{
	string Rep_name = "RepDays=" + std::to_string(Params::Num_Rep_Days) + ".txt";
	vector<int> Tg;  // days of planning
	vector<int> RepDays;
	vector<int> RepDaysCount;
	Read_rep_days(Rep_name, RepDays, RepDaysCount);

	const int nRepDays = (int)RepDays.size();
	int PP = 24 * nRepDays;  // hours of planning for electricity network
	vector<int> Te;
	vector<int> time_weight;
	for (int i = 0; i < nRepDays; i++)
	{
		for (int j = 0; j < 24; j++)
		{
			Te.push_back(24 * RepDays[i] + j);
			time_weight.push_back(RepDaysCount[i]);
			if (Te.size() >= PP) { break; }
		}
	}
	for (int i = 0; i < nRepDays; i++)
	{
		Tg.push_back(RepDays[i]);
	}
	std::map<int, vector<int>> Lg; //key: from_ng_node*200+to_ng_node, 200 is up_lim for number of buses
	vector<gnode> Gnodes = gnode::read_gnode_data("ng_nodes.txt");
	gnode::read_Lexp_data("ng_L_exp.txt", Gnodes);
	gnode::read_Limp_data("ng_L_imp.txt", Gnodes);

	gnode::read_adjE_data("ng_adjE.txt", Gnodes);
	gnode::read_ng_demand_data("ng_daily_dem.txt", Gnodes);



	vector<exist_gSVL> Exist_SVL = exist_gSVL::read_exist_SVL_data("ng_exist_SVL.txt");
	vector<SVL> SVLs = SVL::read_SVL_data("SVL_data.txt");
	vector<pipe> PipeLines = pipe::read_pipe_data("g2g_br.txt", Lg);
	int nGnode = (int)Gnodes.size();

	// Read Electricity Data
	// better to use std::unordered_map which is more efficient and faster
	std::map<int, vector<int>> Le; //key: from_bus*200+to_bus, 200 is up_lim for number of buses
	vector<eStore> Estorage = eStore::read_elec_storage_data("storage_elec.txt");
	vector<enode> Enodes = enode::read_bus_data("bus_num.txt");
	enode::read_adj_data("bus_adj_Nodes.txt", Enodes);
	enode::read_exist_plt_data("existing_plants.txt", Enodes);
	enode::read_demand_data("elec_dem_per_zone_per_hour.txt", Enodes);
	int nEnode = (int)Enodes.size();

	vector<plant> Plants = plant::read_new_plant_data("plant_data.txt");
	plant::read_VRE_profile("profile_hydro_hourly.txt",
		"profile_wind_hourly.txt", "profile_solar_hourly.txt", Plants);

	vector<branch> Branches = branch::read_branch_data(nEnode, "Branches.txt", Le);




	Params::Tg = Tg;
	Params::Te = Te;
	Params::time_weight = time_weight;
	Params::RepDaysCount = RepDaysCount;
	Params::Branches = Branches;
	Params::Enodes = Enodes;
	Params::Gnodes = Gnodes;
	Params::Exist_SVL = Exist_SVL;
	Params::SVLs = SVLs;
	Params::PipeLines = PipeLines;
	Params::Plants = Plants;
	Params::Estorage = Estorage;
	Params::Le = Le;
	Params::Lg = Lg;
}