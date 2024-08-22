#pragma once

#include "vectorOperation.h"
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <random>
#include <chrono>
#include <set>
#include <omp.h>

const double pi = 3.1416;
#define BAND3050

// 2~2.5um band
#ifdef BAND2025
const int GN1 = 15;
const int M1 = 3196;
const int chosen_rows_num = 20;

const int GN2 = 5;
const int M2 = 728;
const int chosen_cols_num = 20;

const double lower_wave_number = 4000.0;
const double upper_wave_number = 5000.0;
#endif

// 3~5um band
#ifdef BAND3050
const int GN1 = 6;
const int M1 = 16;
const int chosen_rows_num = 4;

const int GN2 = 13;
const int M2 = 96;
const int chosen_cols_num = 4;

const double lower_wave_number = 2000.0;
const double upper_wave_number = 3333.0;
#endif

// 3.7~4.8um band
#ifdef BAND3748
const int GN1 = 5;
const int M1 = 247;
const int chosen_rows_num = 20;

const int GN2 = 10;
const int M2 = 10000;
const int chosen_cols_num = 20;

const double lower_wave_number = 2083.0;
const double upper_wave_number = 2703.0;
#endif

// 7.7~9.7um band
#ifdef BAND7797
const int GN1 = 11;
const int M1 = 10000;
const int chosen_rows_num = 40;

const int GN2 = 2;
const int M2 = 10;
const int chosen_cols_num = 10;

const double lower_wave_number = 1031.0;
const double upper_wave_number = 1299.0;
#endif

// 8~14um band
#ifdef BAND0814
const int GN1 = 10;
const int M1 = 10000;
const int chosen_rows_num = 20;

const int GN2 = 10;
const int M2 = 2172;
const int chosen_cols_num = 20;

const double lower_wave_number = 715.0;
const double upper_wave_number = 1250.0;
#endif


const int parallel_threads_number = 32;
const int testN = 56;
const int output_point_number = 401;
const int g_point_number2 = 100000;
const double dyita = 0.005;

const int ngls1 = 3; // H2O高斯积分点个数起始值
const int ngls2 = 3; // CO2高斯积分点个数起始值
const int tng1 = 7; // H2O高斯积分点取值个数
const int tng2 = 7; // CO2高斯积分点取值个数
const int tnp1 = 9; // H2O参考温度取值个数
const int tnp2 = 9; // CO2参考温度取值个数
const int dng = 1; // 高斯积分点取值间隔, both


double*** interval_KT1;
int** interval_KG2;

int interval_number0, interval_number1;
int gn[2], startN[testN][2], NM, NMn[testN];

double TP[tnp1 + tnp2];
double output_parameter2[testN][output_point_number];
int error1[GN1][tng1][tnp1][testN][output_point_number], error2[GN2][tng2][tnp2][testN][output_point_number];
double aerosol_35[4][12001], mole_concentration[6], PM[4][16][6], LH0;


void M5(int j, double pressure)
{
	double b0, b1;
	int i = 0, k;

	while (pressure < PM[j][i][0] && i < 16)
	{
		i++;
	}

	if (i == 0)
		for (k = 1; k < 6; k++)
			mole_concentration[k] = PM[j][0][k];
	else if (i == 16)
		for (k = 1; k < 6; k++)
			mole_concentration[k] = PM[j][15][k];
	else
	{
		b1 = (pressure - PM[j][i][0]) / (PM[j][i - 1][0] - PM[j][i][0]);
		for (k = 1; k < 6; k++)
			mole_concentration[k] = b1 * PM[j][i - 1][k] + (1.0 - b1) * PM[j][i][k];
	}

	LH0 = 0.3 / pressure;
}

// 个体
struct Individual {

	vec1i vars; // 变量
	int error; // 计算误差
	int eqts; // 积分点总数
	double crowding_distance; // 拥挤度

	Individual() {
		this->vars = vec1i{};
		this->error = 100000000000;
		this->eqts = 100000000000;
		this->crowding_distance = -1;
	}

	Individual(int vars_num) {
		this->vars = vec1i(vars_num, -1);
		this->error = 100000000000;
		this->eqts = 100000000000;
		this->crowding_distance = -1;
	}

};


struct VGA {

	int pop_size; // 种群规模
	int max_gen;  // 最大迭代次数
	double cross_pr; // 交叉概率
	double mutation_pr; // 变异概率
	int mutation_num; // 变异点个数

	int Gauss_points_num1; // H2O高斯积分点个数
	int Gauss_points_num2; // CO2高斯积分点个数
	int Gauss_vals_num1; // H2O高斯积分点取法数量
	int Gauss_vals_num2; // CO2高斯积分点取法数量
	int Tp_num1;  // H2O的参考温度个数
	int Tp_num2;  // CO2的参考温度个数
	int total_pos_num; // 总变量有多少位

	int strategy_idx1; // H2O的方案编号
	int strategy_idx2; // H2O的方案编号

	vec1i val_range_map; // 存储对应索引位置的变量类型，1：H2O高斯积分点数，2：CO2高斯积分点数，3：H2O参考温度，4：CO2参考温度

	int father, son; // 指示当前种群参数存在0还是1

	vector<vector<Individual> > individuals; // 种群所有个体

	VGA(int pop_size, int max_gen, int Gauss_points_num1, int Gauss_points_num2, int Gauss_vals_num1, int Gauss_vals_num2, int strategy_idx1, int strategy_idx2, int Tp_num1, int Tp_num2) {

		this->pop_size = pop_size;
		this->max_gen = max_gen;
		this->cross_pr = 0.8;
		this->mutation_pr = 0;
		this->mutation_num = 1;

		this->Gauss_points_num1 = Gauss_points_num1;
		this->Gauss_points_num2 = Gauss_points_num2;
		this->Tp_num1 = Tp_num1;
		this->Tp_num2 = Tp_num2;
		this->Gauss_vals_num1 = Gauss_vals_num1;
		this->Gauss_vals_num2 = Gauss_vals_num2;
		this->strategy_idx1 = strategy_idx1;
		this->strategy_idx2 = strategy_idx2;

		this->total_pos_num = 2 * (Gauss_points_num1 + Gauss_points_num2); // 总变量有多少位：H2O高斯积分点数+H2O参考温度+CO2高斯积分点数+CO2参考温度+H2O分组方案

		// 初始化val_range_map
		this->val_range_map = vec1i(total_pos_num, 0);
		// H2O
		for (int i = 0; i < Gauss_points_num1; ++i) {
			val_range_map[2 * i] = 1; // H2O高斯积分点数
			val_range_map[2 * i + 1] = 3; // H2O参考温度
		}
		// CO2
		for (int i = 0; i < Gauss_points_num2; ++i) {
			val_range_map[2 * Gauss_points_num1 + 2 * i] = 2; // CO2高斯积分点数
			val_range_map[2 * Gauss_points_num1 + 2 * i + 1] = 4; // CO2参考温度
		}

		this->father = 0;
		this->son = 1;

		this->individuals = vector<vector<Individual> >(2, vector<Individual>(pop_size, Individual(total_pos_num)));
	}
};

// 根据不同基因点位生成合适的随机数
inline int rand_val(VGA& data, int index) {

	const vec1i& val_range_map = data.val_range_map;

	// 生成整数的随机值引擎
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr_int_gauss1(0, data.Gauss_vals_num1 - 1); // H2O高斯积分点的索引，随机整数
	uniform_int_distribution<int> distr_int_gauss2(0, data.Gauss_vals_num2 - 1); // CO2高斯积分点的索引，随机整数
	uniform_int_distribution<int> distr_int_Tp1(0, data.Tp_num1 - 1); // Tp1的索引，随机整数
	uniform_int_distribution<int> distr_int_Tp2(0, data.Tp_num2 - 1); // Tp2的索引，随机整数

	if (val_range_map[index] == 1)
		return distr_int_gauss1(eng);
	else if (val_range_map[index] == 2)
		return distr_int_gauss2(eng);
	else if (val_range_map[index] == 3)
		return distr_int_Tp1(eng);
	else if (val_range_map[index] == 4)
		return distr_int_Tp2(eng);
	else {
		printf("[ERROR] rand_val: index = %d...!n", index);
		exit(-1);
	}

}

int init_GA(VGA& data) {

	const int& father = data.father;

	// 第一代父种群变量赋随机值
	for (int i = 0; i < data.pop_size; ++i) {
		for (int j = 0; j < data.total_pos_num; ++j)
			data.individuals[father][i].vars[j] = rand_val(data, j);
	}

	return 0;
}


bool cmp_by_error(Individual& a, Individual& b) {
	return a.error < b.error;
}

bool cmp_by_eqts(Individual& a, Individual& b) {
	return a.eqts < b.eqts;
}

bool cmp_by_crowding_distance(Individual& a, Individual& b) {
	return a.crowding_distance > b.crowding_distance; // 注意是大于
}

// 拥挤度排序用的结构体
struct Front_Points {

	int index; // 在data中的索引值
	int error; // 计算误差
	int eqts; // 积分点总数
	double crowding_distance;

	Front_Points(Individual& individual, int i) {
		this->index = i;
		this->error = individual.error;
		this->eqts = individual.eqts;
		this->crowding_distance = 0;
	}

};

bool cmp_by_error1(Front_Points& a, Front_Points& b) {
	return a.error < b.error;
}

bool cmp_by_eqts1(Front_Points& a, Front_Points& b) {
	return a.eqts < b.eqts;
}

bool cmp_by_crowding_distance1(Front_Points& a, Front_Points& b) {
	return a.crowding_distance > b.crowding_distance; // 注意是大于
}

vec1i picked_by_crowding_distance(VGA& data, vec1i& front, int num) {

	vec1i picked_nums(num);

	int& father = data.father;
	int& son = data.son;
	double delta;

	vector<Front_Points> points{}; // 当前front的所有个体
	for (int i = 0; i < front.size(); ++i) {
		//points.emplace_back(front[i], data.individuals[father][front[i]].error, data.individuals[father][front[i]].eqts, 0);
		points.emplace_back(data.individuals[father][front[i]], front[i]);
	}


	// 按error从小到大排序
	sort(points.begin(), points.end(), cmp_by_error1);
	delta = points[points.size() - 1].error - points[0].error; // 最大error和最小error的差
	if (delta != 0) { // 分母不为0

		// 边缘两个拥挤度无限大
		points[0].crowding_distance += 1e20;
		points[points.size() - 1].crowding_distance += 1e20;

		for (int i = 1; i < points.size() - 1; ++i)
			points[i].crowding_distance += (points[i + 1].error - points[i - 1].error) / delta;
	}

	// 按eqts从小到大排序
	sort(points.begin(), points.end(), cmp_by_eqts1);
	delta = points[points.size() - 1].eqts - points[0].eqts; // 最大eqts和最小eqts的差
	if (delta != 0) { // 分母不为0

		// 边缘两个拥挤度无限大
		points[0].crowding_distance += 1e20;
		points[points.size() - 1].crowding_distance += 1e20;

		for (int i = 1; i < points.size() - 1; ++i)
			points[i].crowding_distance += (points[i + 1].eqts - points[i - 1].eqts) / delta;
	}

	// 按算出的拥挤度从大到小排序
	sort(points.begin(), points.end(), cmp_by_crowding_distance1);


	// 按照剩余名额挑出拥挤度最大的个体
	for (int i = 0; i < num; ++i) {
		picked_nums[i] = points[i].index;
	}

	return picked_nums;
}

// p支配q为真
inline bool is_dominating(VGA& data, int p, int q) {

	int& p1 = data.individuals[data.father][p].eqts;
	int& q1 = data.individuals[data.father][q].eqts;
	int& p2 = data.individuals[data.father][p].error;
	int& q2 = data.individuals[data.father][q].error;

	if ((p1 < q1 && p2 < q2) ||
		(p1 == q1 && p2 < q2) ||
		(p1 < q1 && p2 == q2))
		return true;
	else
		return false;
}

// 计算适应度
int calc_error(VGA& data, int beg_ind, int end_ind, int startN[][2], const vec1i& gauss_map1, const vec1i& gauss_map2) {

	int b0, b1, b2, b3, b4, x, n3, i, j, k;
	int& father = data.father;
	int error00[testN];
#pragma omp parallel for schedule(static) private(i,j,k,x,n3,b0,b1,b2,b3,b4,error00)
	for (n3 = 0; n3 < parallel_threads_number; n3++)
	{
		for (x = beg_ind + n3; x < end_ind; x += parallel_threads_number)
		{
			Individual& individual = data.individuals[father][x];
			vec1i& vars = individual.vars;
			i = testN - 1;
			do {

				for (b0 = 0, j = startN[i][0]; j <= startN[i][1]; ++j)
				{
					b1 = 0;
					for (k = 0; k < GN2; k++) {
						//printf("%d, %d, %d, %d, %d, %d\n", n3 / threads_per_numa, k, vars[k + GN1], vars[data.Tp_index2], i, j);
						b1 += error2[k][vars[2 * GN1 + 2 * k]][vars[2 * GN1 + 2 * k + 1]][i][j];
					}

					for (k = 0; k < GN1; k++) {
						b1 += error1[k][vars[2 * k]][vars[2 * k + 1]][i][j];
					}
					if (abs(b1) > abs(b0))
						b0 = b1;
				}
				error00[i] = b0;
				i--;
			} while (i >= 0);

			for (b4 = 0, i = 0; i < testN; i++)
			{
				b0 = abs(error00[i]);

				b2 = b0 - 800;
				b2 = b2 > 0 ? b2 * 9 : 0;

				b3 = b0 - 1200;
				b3 = b3 > 0 ? b3 * 90 : 0;

				b4 += (b0 + b2 + b3);
			}
			// 计算两个目标函数值
			individual.error = b4;  // 该个体的误差值
			individual.eqts = 0;
			for (int ii = 0; ii < data.Gauss_points_num1; ++ii) // 注意每次自增2
				individual.eqts += gauss_map1[individual.vars[2 * ii]];
			for (int ii = 0; ii < data.Gauss_points_num2; ++ii) // 注意每次自增2
				individual.eqts += gauss_map2[individual.vars[2 * data.Gauss_points_num1 + 2 * ii]];
		}
	}

	return 0;
}

// 给出需要计算的个体的下标集合vec_index，只计算这些个体的适应度
int calc_error2(VGA& data, vec1i vec_index, int startN[][2], const vec1i& gauss_map1, const vec1i& gauss_map2) {

	int b0, b1, b2, b3, b4, x, n3, i, j, k;
	int& father = data.father;
	int error00[testN];
#pragma omp parallel for schedule(static) private(i,j,k,x,n3,b0,b1,b2,b3,b4,error00)
	for (n3 = 0; n3 < parallel_threads_number; n3++)
	{
		for (int ii = n3; ii < vec_index.size(); ii += parallel_threads_number) {

			x = vec_index[ii];
			Individual& individual = data.individuals[father][x];
			vec1i& vars = individual.vars;
			i = testN - 1;
			do {

				for (b0 = 0, j = startN[i][0]; j <= startN[i][1]; ++j)
				{
					b1 = 0;
					for (k = 0; k < GN2; k++) {
						b1 += error2[k][vars[2 * GN1 + 2 * k]][vars[2 * GN1 + 2 * k + 1]][i][j];
					}

					for (k = 0; k < GN1; k++) {
						b1 += error1[k][vars[2 * k]][vars[2 * k + 1]][i][j];
					}
					if (abs(b1) > abs(b0))
						b0 = b1;
				}
				error00[i] = b0;
				i--;
			} while (i >= 0);

			for (b4 = 0, i = 0; i < testN; i++)
			{
				b0 = abs(error00[i]);

				b2 = b0 - 800;
				b2 = b2 > 0 ? b2 * 9 : 0;

				b3 = b0 - 1200;
				b3 = b3 > 0 ? b3 * 90 : 0;

				b4 += (b0 + b2 + b3); // 这里没有乘ww了
			}
			// 计算两个目标函数值
			individual.error = b4;  // 该个体的误差值
			individual.eqts = 0;
			for (int ii = 0; ii < data.Gauss_points_num1; ++ii) // 注意每次自增2
			{
				individual.eqts += gauss_map1[individual.vars[2 * ii]];

			}
			for (int ii = 0; ii < data.Gauss_points_num2; ++ii) // 注意每次自增2
			{
				individual.eqts += gauss_map2[individual.vars[2 * data.Gauss_points_num1 + 2 * ii]];
			}
		}
	}

	return 0;
}

// generate_offspring_NSGA2专用
int non_dominated_sort_2(VGA& data, int chosen_num, int startN[][2], const vec1i& gauss_map1, const vec1i& gauss_map2, int max_mutation_repeat_times, int repeat_times) {

	// 生成整数的随机值引擎
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr_int_mutation(0, data.total_pos_num - 1); // 变异选点
	uniform_int_distribution<int> distr_int_mutation_gauss(0, data.Gauss_points_num1 + data.Gauss_points_num2 - 1); // 变异选点
	uniform_real_distribution<double> distr_double(0, 1); // 0~1概率，随机浮点数


	int& father = data.father;
	//int& H2O_index = data.H2O_index;
	vector<Individual>& individuals = data.individuals[father];

	sort(individuals.begin(), individuals.end(), cmp_by_error);

	//if (repeat_times > 5) {
	int err_now = individuals[0].error;
	int eqts_now = individuals[0].eqts;

	int ind_mutation = -1; // 变异在个体的哪个变量位
	int mutation_bit_num = -1; // 变异点位数量
	vec1i vec_index{}; // 保存需要重新计算目标函数的个体的索引
	for (int i = 1; i < individuals.size(); ++i) {
		if (err_now != individuals[i].error || eqts_now != individuals[i].eqts) {
			err_now = individuals[i].error;
			eqts_now = individuals[i].eqts;
		}
		else { // 随机取两位做变异
			//if (distr_double(eng) < 0.5) {
				//mutation_bit_num = distr_double(eng) < 0.5 ? 1 : 2; // 以均等概率变异1或2个点
				//for (int j = 0; j < mutation_bit_num; ++j) {
				//	ind_mutation = distr_int_mutation(eng);
				//	//ind_mutation = data.Tp_index_vec[ind_mutation];
				//	individuals[i].vars[ind_mutation] = rand_val(data, ind_mutation);
				//}
			//	ind_mutation = distr_int_mutation(eng);
			//	individuals[i].vars[ind_mutation] = rand_val(data, ind_mutation);
			//}
			//else {
			ind_mutation = distr_int_mutation_gauss(eng);
			individuals[i].vars[2 * ind_mutation] = rand_val(data, 2 * ind_mutation);
			individuals[i].vars[2 * ind_mutation + 1] = rand_val(data, 2 * ind_mutation + 1);
			//}
			vec_index.push_back(i);

		}
	}
	calc_error2(data, vec_index, startN, gauss_map1, gauss_map2);
	//}


	vec2i fronts(data.pop_size, vec1i{}); // 帕雷托前沿
	vec1i rank(data.pop_size, 0); // 等级
	vec2i S(data.pop_size, vec1i{}); // p的支配集合
	vec1i N(data.pop_size, 0); // p的被支配度

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < data.pop_size; ++i) {
		for (int j = 0; j < data.pop_size; ++j) {
			if (is_dominating(data, i, j)) { // 如果i支配j
				S[i].push_back(j); // j纳入i的支配集合
			}
			else if (is_dominating(data, j, i)) { //如果i被j支配
				N[i] += 1; //i的被支配度增1
			}
		}
		if (N[i] == 0)
			rank[i] = 0;
	}

	for (int i = 0; i < data.pop_size; ++i) {
		if (N[i] == 0)
			fronts[0].push_back(i);
	}

	int total_num = fronts[0].size(); // 被选中的前沿点的总数
	int k = 0;
	vec1i Q{};
	while (fronts[k].size() != 0 && total_num < chosen_num) {
		Q.clear();
		for (int i = 0; i < fronts[k].size(); ++i) { // 当前最前沿的每个点
			int& idx = fronts[k][i]; // 当前前沿点的实际编号
			//cout << "------" << idx << endl;
			for (int j = 0; j < S[idx].size(); ++j) { // 每个点支配的点的集合
				int& ind = S[idx][j]; // 被支配的点的实际索引号
				//cout << ind << "/ " << endl;
				N[ind] -= 1; // 当前前沿点去除后
				if (N[ind] == 0) { // ，该被支配点能否成为新的前沿点？
					rank[ind] = k + 1; // rank + 1
					Q.push_back(ind); // 该被支配点被加入下一组前沿点
				}
			}
		}
		++k;
		fronts[k] = Q;
		total_num += Q.size(); // 记录目前已经选中的前沿点的总数
	}

	int& son = data.son;

	// 挑选按拥挤度生成新的父代
	int cnt = 0; // 目前新一代个体数
	for (int i = 0; i < k + 1; ++i) {
		if (cnt == chosen_num)
			break;

		// 目前新一代个体数加上新前沿个体数仍不大于pop_size
		if (cnt + fronts[i].size() <= chosen_num) {
			for (int j = 0; j < fronts[i].size(); ++j) {
				int& ind = fronts[i][j];
				data.individuals[son][cnt] = data.individuals[father][ind];
				++cnt;
			}
		}
		// 目前新一代个体数加上新前沿个体数大于pop_size，计算拥挤度挑选
		else {
			int left_num = chosen_num - cnt; // 剩余名额
			vec1i picked_nums = picked_by_crowding_distance(data, fronts[i], left_num);
			for (int j = 0; j < left_num; ++j) {
				int& ind = picked_nums[j];
				data.individuals[son][cnt] = data.individuals[father][ind];
				++cnt;
			}
		}
	}


	// 检查cnt是不是等于pop_size
	if (cnt != chosen_num) {
		int left_num = chosen_num - cnt;
		cout << "cnt == " << cnt << ", chosen_num = " << chosen_num << endl;
		exit(-1);
	}

	return 0;
}

int generate_offspring_NSGA2(VGA& data, int GA_sgn, int beg_ind, int end_ind, int startN[][2], const vec1i& gauss_map1, const vec1i& gauss_map2, int max_mutation_repeat_times, int repeat_times) {
	double& Pc = data.cross_pr; // 交叉概率
	double& Pm = data.mutation_pr; // 变异概率
	//int& H2O_index = data.H2O_index;
	const int& total_pos_num = data.total_pos_num;

	random_device rd;
	default_random_engine eng(rd());
	uniform_real_distribution<double> distr_double(0, 1); // 0~1概率，随机浮点数
	uniform_int_distribution<int> distr_int_pop(0, data.pop_size - 1); // 种群数，随机整数
	uniform_int_distribution<int> distr_int_cross(0, total_pos_num); // 交叉选点
	uniform_int_distribution<int> distr_int_mutation_num(1, data.mutation_num); // 变异选点-total
	uniform_int_distribution<int> distr_int_mutation(0, total_pos_num - 1); // 变异选点-total


	int& father = data.father; // 当前父辈
	int& son = data.son;

	int ind1, ind2; // 两个被选择的父代编号
	double tmp1, tmp2;
	int ind_cross1, ind_cross2, ind_cross3, ind_mutation;
	double roll_val = -1;
	int mutation_point_num = -1; // 确定变异点有几个

	int half_size = data.pop_size / 2;

	calc_error(data, half_size, data.pop_size, startN, gauss_map1, gauss_map2);

	// 非支配排序，取前half_size进入子代
	non_dominated_sort_2(data, half_size, startN, gauss_map1, gauss_map2, max_mutation_repeat_times, repeat_times);

	// 对子代的half_size个个体交叉变异生成另一半子代
		// Pr_table：轮盘赌用的数组, pop_size+1项, 0到1
	double err_max = -1e10;
	vec1d Pr_table(half_size + 1, 0);
	for (int i = 0; i < half_size; ++i) {
		err_max = data.individuals[son][i].error > err_max ? data.individuals[son][i].error : err_max;
		Pr_table[i + 1] = data.individuals[son][i].error + Pr_table[i]; // 0815更新，之前是-data.errors[i] + Pr_table[i]，删去了负号
	}

	// 如果所有个体都一样了，已经收敛，不用再算了
	if (half_size * err_max - *(Pr_table.end() - 1) == 0) {
		return 1;
	}

	for (int i = 0; i < half_size; ++i) {
		Pr_table[i] = (i * err_max - Pr_table[i]) / (half_size * err_max - *(Pr_table.end() - 1));
	}

	// 考虑最小化eqts，加权一个轮盘赌
	double eqts_max = -1;
	vec1d Eq_table(half_size + 1, 0);
	for (int i = 0; i < half_size; ++i) {
		err_max = data.individuals[son][i].eqts > err_max ? data.individuals[son][i].eqts : err_max;
		Eq_table[i + 1] = data.individuals[son][i].eqts + Eq_table[i]; // 0815更新，之前是-data.errors[i] + Pr_table[i]，删去了负号
	}
	for (int i = 0; i < half_size; ++i) {
		Eq_table[i] = (i * eqts_max - Eq_table[i]) / (half_size * eqts_max - *(Eq_table.end() - 1));
		Pr_table[i] = (Pr_table[i] + Eq_table[i]) / 2;
	}

	//vec1i break_points(3);

	// 交叉和变异
	for (int i = half_size; i < data.pop_size; ++i) { // 从elite_num之后，产生新的子代种群

		tmp1 = distr_double(eng);
		tmp2 = distr_double(eng);
		if (tmp1 == 0)
			ind1 = 0;
		else if (tmp1 == *(Pr_table.end() - 1))
			ind1 = data.pop_size - 1;
		else
			for (int j = 0; j < Pr_table.size(); ++j) {
				if (tmp1 < Pr_table[j]) {
					ind1 = j - 1;
					break;
				}
			}
		if (tmp2 == 0)
			ind2 = 0;
		else if (tmp2 == *(Pr_table.end() - 1))
			ind2 = data.pop_size - 1;
		else
			for (int j = 0; j < Pr_table.size(); ++j) {
				if (tmp2 < Pr_table[j]) {
					ind2 = j - 1;
					break;
				}
			}


		if (distr_double(eng) < Pc) { // 两点交叉
			ind_cross1 = distr_int_cross(eng);
			ind_cross2 = distr_int_cross(eng);
			if (ind_cross1 > ind_cross2) {
				int tmp = ind_cross1;
				ind_cross1 = ind_cross2;
				ind_cross2 = tmp;
			}
			for (int x = 0; x < ind_cross1; ++x) {
				data.individuals[son][i].vars[x] = data.individuals[son][ind1].vars[x];
			}
			for (int x = ind_cross1; x < ind_cross2; ++x) {
				data.individuals[son][i].vars[x] = data.individuals[son][ind2].vars[x];
			}
			for (int x = ind_cross2; x < total_pos_num; ++x) {
				data.individuals[son][i].vars[x] = data.individuals[son][ind1].vars[x];
			}
		}
		else { // 直接复制
			for (int x = 0; x < total_pos_num; ++x)
				data.individuals[son][i].vars[x] = data.individuals[son][ind1].vars[x];
		}

		// 计算有几个变异点
		mutation_point_num = distr_int_mutation_num(eng);
		// 全体变异
		for (int ii = 0; ii < mutation_point_num; ++ii) {
			if (distr_double(eng) < Pm) {
				ind_mutation = distr_int_mutation(eng);
				data.individuals[son][i].vars[ind_mutation] = rand_val(data, ind_mutation);
			}
		}
	}
	return GA_sgn;
}


inline int convert_gene_val(const VGA& data, const Individual& individual, const vec1i& gauss_map1, const vec1i& gauss_map2, vec1i& dna) {

	dna = vec1i(individual.vars.size());

	int sgn = -1;
	for (int j = 0; j < individual.vars.size(); ++j) {
		sgn = data.val_range_map[j];
		if (sgn == 1) { // H2O Gauss
			dna[j] = gauss_map1[individual.vars[j]];
		}
		else if (sgn == 2) { // CO2 Gauss
			dna[j] = gauss_map2[individual.vars[j]];
		}
		else
			dna[j] = individual.vars[j];
	}
}

struct Goal {

	vec1i vars; // 所有变量，包括H2O
	double error; // 计算误差
	int eqts; // 积分点总数
	int total_pos_num;
	int H2O_group_num;

	//Goal(int num) {
	//	this->vars = vec1i(num, -1);
	//	this->error = 1000000000;
	//	this->eqts = 1000000000;
	//	this->H2O_group_index = -1;
	//}

	Goal() {
		this->vars = vec1i{};
		this->error = 1000000000;
		this->eqts = 1000000000;
		this->H2O_group_num = -1;
	}

	Goal(VGA& data, Individual& individual, vec1i& gauss_map1, vec1i& gauss_map2) {
		this->vars = vec1i(data.total_pos_num);
		convert_gene_val(data, individual, gauss_map1, gauss_map2, this->vars);
		this->error = individual.error / 10000.0;
		this->eqts = individual.eqts;
		this->H2O_group_num = data.strategy_idx1;
		this->total_pos_num = data.total_pos_num;
	}

};

// 挑选当代最优个体
Goal find_goal(VGA& data, vec1i& gauss_map1, vec1i& gauss_map2) {

	const int& father = data.father;
	int ind = -1;
	double err = 1e20;

	for (int i = 0; i < data.pop_size; ++i) {
		if (data.individuals[father][i].error < err) {
			err = data.individuals[father][i].error;
			ind = i;
		}
	}

	Individual& best = data.individuals[father][ind];

	return Goal(data, best, gauss_map1, gauss_map2);
}

// 交换父子代
int swap_father_son_sgn(VGA& data) {
	int temp = data.father;
	data.father = data.son;
	data.son = temp;
	return 0;
}

struct Pareto_Fronts {
	int solution_num = -1;
	vec2d errors;
	vec2i eqts;

	Pareto_Fronts(int num) {
		this->solution_num = num;
		this->errors = vec2d{};
		this->eqts = vec2i{};

	}
};

struct mat0 : Individual {
	int index;
	mat0(Individual& individual, int i) {
		this->index = i;
		this->vars = individual.vars;
		this->error = individual.error;
		this->eqts = individual.eqts;
	}
};

inline int remove_repeat_point(VGA& data, vec1i& front_index) {

	vector<mat0> front{};
	for (int i = 0; i < front_index.size(); ++i)
		front.emplace_back(data.individuals[data.father][front_index[i]], front_index[i]);

	if (front.size() < 2)
		return 0;

	sort(front.begin(), front.end(), cmp_by_error);

	int err_now = front[0].error;
	int eqts_now = front[0].eqts;

	front_index = { front[0].index };

	for (int i = 1; i < front.size(); ++i)
		if (err_now != front[i].error || eqts_now != front[i].eqts) {
			front_index.push_back(front[i].index);
			err_now = front[i].error;
			eqts_now = front[i].eqts;
		}

	return 0;
}

vector<vector<Goal> > get_fronts(VGA& data, int solution_num, vec1i& gauss_map1, vec1i& gauss_map2) {

	vec2i fronts(data.pop_size, vec1i{}); // 帕雷托前沿
	vec1i rank(data.pop_size, 0); // 等级
	vec2i S(data.pop_size, vec1i{}); // p的支配集合
	vec1i N(data.pop_size, 0); // p的被支配度

	int& chosen_num = solution_num;

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < data.pop_size; ++i) {
		for (int j = 0; j < data.pop_size; ++j) {
			if (is_dominating(data, i, j)) { // 如果i支配j
				S[i].push_back(j); // j纳入i的支配集合
			}
			else if (is_dominating(data, j, i)) { //如果i被j支配
				N[i] += 1; //i的被支配度增1
			}

		}
		if (N[i] == 0)
#pragma omp critical
		{
			rank[i] = 0;
			fronts[0].push_back(i);
		}
	}

	int k = 0;
	vec1i Q{};
	while (fronts[k].size() != 0) {
		Q.clear();
		for (int i = 0; i < fronts[k].size(); ++i) { // 当前最前沿的每个点
			int& idx = fronts[k][i]; // 当前前沿点的实际编号
			for (int j = 0; j < S[idx].size(); ++j) { // 每个点支配的点的集合
				int& ind = S[idx][j]; // 被支配的点的实际索引号
				N[ind] -= 1; // 当前前沿点去除后
				if (N[ind] == 0) { // ，该被支配点能否成为新的前沿点？
					rank[ind] = k + 1; // rank + 1
					Q.push_back(ind); // 该被支配点被加入下一组前沿点
				}
			}
		}
		++k;
		fronts[k] = Q;
	}

	vector<vector<Goal> > res{};

	// 挑选按拥挤度生成新的父代
	int cnt = 0; // 目前新一代个体数
	vec1d errors_tmp;
	vec1i eqts_tmp;
	for (int i = 0; i < k + 1; ++i) {

		if (cnt == chosen_num)
			break;

		vector<Goal> individuals_tmp{};

		remove_repeat_point(data, fronts[i]);

		// 目前新一代个体数加上新前沿个体数仍不大于pop_size
		if (cnt + fronts[i].size() <= chosen_num) {
			for (int j = 0; j < fronts[i].size(); ++j) {
				int& ind = fronts[i][j];
				individuals_tmp.push_back(Goal(data, data.individuals[data.father][ind], gauss_map1, gauss_map2));
				++cnt;
			}
		}
		// 目前新一代个体数加上新前沿个体数大于pop_size，计算拥挤度挑选
		else {
			int left_num = chosen_num - cnt; // 剩余名额
			vec1i picked_nums = picked_by_crowding_distance(data, fronts[i], left_num);
			for (int j = 0; j < left_num; ++j) {
				int& ind = picked_nums[j];
				individuals_tmp.push_back(Goal(data, data.individuals[data.father][ind], gauss_map1, gauss_map2));
				++cnt;
			}
		}

		res.push_back(individuals_tmp);

	}

	// 检查cnt是不是等于pop_size
	//if (cnt != chosen_num) {
	//	cout << "[ERROR] cnt == " << cnt << ", chosen_num = " << chosen_num << endl;
	//	exit(-1);
	//}

	return res;
}


int NSGA2_one(int strategy_idx1, int strategy_idx2, int pop_size, int max_gen, string method_name, int startN[][2], Goal& err_now) {
	// 遗传算法begin

	double total_it = 0; // 记录总迭代次数
	double total_time = 0; // 记录总时间

	// 计时开始
	auto start = chrono::steady_clock::now();

	double error_min = 1e20;

	const int vars_num = GN1 + GN2; // 变量数量：水蒸气GN1位+二氧化碳GN2位+水蒸气分组方案1位
	const int max_muation_bit = int(vars_num / 2); // 最大变异位数
	constexpr double cross_pr = 0.8; // 单点交叉概率


	// 给定超参数
	VGA data(pop_size, max_gen, GN1, GN2, tng1, tng2, strategy_idx1, strategy_idx2, tnp1, tnp2);

	int elite_num = int(0.05 * pop_size); // 精英个体数量

	// 初始化种群
	init_GA(data);

	vec1i gauss_map1{}, gauss_map2{};
	for (int i = 0; i < tng1; ++i)
		gauss_map1.push_back(ngls1 + dng * i);
	for (int i = 0; i < tng2; ++i)
		gauss_map2.push_back(ngls2 + dng * i);

	calc_error(data, 0, pop_size, startN, gauss_map1, gauss_map2);

	Goal err_best{}; // 刚初始化，现在不可使用

	int repeat_times = 0;

	char folder_name[256];
	snprintf(folder_name, sizeof(folder_name) - 1, ".\\%s", method_name);
	mkdir(folder_name);

	char folder_name_detail[256];
	snprintf(folder_name_detail, sizeof(folder_name_detail) - 1, "%s\\detail-%d", folder_name, pop_size);
	mkdir(folder_name_detail);

	FILE* fp;
	int removed_it = 50; // 去除重复个体循环次数
	int GA_sgn = 3; // 1：GA， 2：NSGA2

	int max_repeat_times = 50; // 达到最优解停滞迭代次数后结束循环
	int max_mutation_repeat_times = int(0.2 * max_repeat_times);
	int max_mutation_repeat_times_2 = int(0.1 * max_repeat_times); // 变异概率随最优解停滞迭代次数增长函数的分母
	int max_mutation_repeat_times_1 = int(0.5 * max_repeat_times); // 变异概率随最优解停滞迭代次数增长函数的分母
	int change_method_times = int(0.1 * max_repeat_times); // 每迭代几次换GA方法
	int only_NSGA2_times = int(0.6 * max_repeat_times); // 超出该迭代次数后只使用NSGA2

	// 遗传算法迭代begin
	for (int it = 0; it < max_gen; ++it, ++repeat_times) {

		printf("%s generation %d, %d, %.3f...\r", method_name, it, repeat_times, data.mutation_pr);


		// 生成子代：计算适应度、选择、交叉、变异
		generate_offspring_NSGA2(data, GA_sgn, 0, pop_size, startN, gauss_map1, gauss_map2, max_mutation_repeat_times, repeat_times);

		// 当前代最优解
		err_best = find_goal(data, gauss_map1, gauss_map2);

		// 如果有更优解
		if (err_best.error < error_min) {

			// 控制台输出
			printf("\n%d || %d || %.2f %d || %d %d || ", strategy_idx2, it, err_best.error, err_best.eqts, err_best.H2O_group_num, strategy_idx2);
			for (int i = 0; i < err_best.total_pos_num; ++i) {
				printf("%d ", err_best.vars[i]);
			}
			printf("\n");
			error_min = err_best.error;

			// 文件输出
			char name[256];
			snprintf(name, sizeof(name) - 1, "%s\\res-%d-%d.txt", folder_name_detail, strategy_idx1, strategy_idx2);
			fp = fopen(name, "a");
			fprintf(fp, "%d %.4f %d %d %d || ", it, err_best.error, err_best.eqts, err_best.H2O_group_num, strategy_idx2);
			for (int i = 0; i < err_best.total_pos_num; ++i) {
				fprintf(fp, "%d ", err_best.vars[i]);
			}
			fprintf(fp, "\n");
			fclose(fp);

			// 最优解停滞的代数归零
			repeat_times = 0;
			//repeat_times = int(0.1 * max_mutation_repeat_times);
		}

		// 如果最优解停滞的代数超过500代，或者循环到了最后一次，写文件并结束迭代
		if (repeat_times > max_repeat_times || it == max_gen - 1) {

			// 文件输出
			char name[256];
			snprintf(name, sizeof(name) - 1, "%s\\res-%d-%d.txt", folder_name_detail, strategy_idx1, strategy_idx2);
			fp = fopen(name, "a");
			fprintf(fp, "%d %.4f %d %d %d || ", it, err_best.error, err_best.eqts, err_best.H2O_group_num, strategy_idx2);
			for (int i = 0; i < err_best.total_pos_num; ++i) {
				fprintf(fp, "%d ", err_best.vars[i]);
			}
			fprintf(fp, "\n");
			fclose(fp);

			total_it = it; // 记录总迭代次数

			break; // 结束迭代
		}


		// 如果最优解停滞的代数超过100代，增加变异点位数量，以增强解的多样性，扩大搜索范围
		else if (repeat_times > 20) {
			data.mutation_num = 2;
			data.mutation_pr = 0.2;
		}
		// 随着最优解停滞的代数的增多，变异概率逐渐增大到1
		else {
			data.mutation_num = 1;
			data.mutation_pr = 0.05;
		}
		// 翻转父代子代标识
		swap_father_son_sgn(data);

	}

	// 遗传算法迭代end
	printf("\n");

	// 写当前组的最优解
	{

		// 写当前最优解2
		{
			char name[256];
			snprintf(name, sizeof(name) - 1, "%s\\object_error-%d.txt", folder_name, pop_size);
			fp = fopen(name, "a");
			fprintf(fp, "%d || %.4f %d || %d %d || ", strategy_idx2, err_best.error, err_best.eqts, err_best.H2O_group_num, strategy_idx2);
			for (int i = 0; i < err_best.total_pos_num; ++i) {
				fprintf(fp, "%d ", err_best.vars[i]);
			}
			fprintf(fp, "\n");
			fclose(fp);
		}

		// 写pareto前沿结果，每个zone代表一层前沿
		{
			vector<vector<Goal> > pareto_fronts = get_fronts(data, elite_num, gauss_map1, gauss_map2);
			char name[256];
			snprintf(name, sizeof(name) - 1, "%s\\fronts-%d-%d.txt", folder_name_detail, strategy_idx1, strategy_idx2);
			fp = fopen(name, "w");
			fprintf(fp, "VARIABLES = \"Quadrature points total quantity\", \"Error value\"\n");

			for (int it_f = 0; it_f < pareto_fronts.size(); ++it_f) {
				fprintf(fp, "ZONE T = \"Front %d\", I = %d, DATAPACKING = POINT\n", it_f + 1, pareto_fronts[it_f].size());
				for (int it_p = 0; it_p < pareto_fronts[it_f].size(); ++it_p) {
					fprintf(fp, "%d %.4f\n", pareto_fronts[it_f][it_p].eqts, pareto_fronts[it_f][it_p].error);
				}
			}
			fprintf(fp, "\n");
			fclose(fp);
		}

		// 写pareto前沿详细结果，每个zone代表一层前沿
		{
			vector<vector<Goal> > pareto_fronts = get_fronts(data, elite_num, gauss_map1, gauss_map2);
			char name[256];
			snprintf(name, sizeof(name) - 1, "%s\\fronts-details-%d-%d.txt", folder_name_detail, strategy_idx1, strategy_idx2);
			fp = fopen(name, "w");
			for (int it_f = 0; it_f < pareto_fronts.size(); ++it_f) {
				fprintf(fp, "Front %d, I = %d\n", it_f + 1, pareto_fronts[it_f].size());
				for (int it_p = 0; it_p < pareto_fronts[it_f].size(); ++it_p) {
					fprintf(fp, "%d %.4f || ", pareto_fronts[it_f][it_p].eqts, pareto_fronts[it_f][it_p].error);
					fprintf(fp, "%d %d || ", pareto_fronts[it_f][it_p].H2O_group_num, strategy_idx2);
					for (int it_v = 0; it_v < data.total_pos_num; ++it_v) {
						fprintf(fp, "%d ", pareto_fronts[it_f][it_p].vars[it_v]);
					}
					fprintf(fp, "\n");
				}
			}
			fprintf(fp, "\n");
			fclose(fp);
		}
	}

	// 计时结束;
	auto end = chrono::steady_clock::now();
	auto time_diff = end - start;
	auto duration = chrono::duration_cast<chrono::seconds>(time_diff);
	cout << "Operation cost : " << duration.count() << "s" << endl;

	total_time = chrono::duration_cast<chrono::milliseconds>(time_diff).count() / 1e3; // 记录总时间(s)

	err_now = err_best;

	return 0;
	//return vec1d{ total_it, total_time, total_time / total_it };
}

