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

const int ngls1 = 3; // H2O��˹���ֵ������ʼֵ
const int ngls2 = 3; // CO2��˹���ֵ������ʼֵ
const int tng1 = 7; // H2O��˹���ֵ�ȡֵ����
const int tng2 = 7; // CO2��˹���ֵ�ȡֵ����
const int tnp1 = 9; // H2O�ο��¶�ȡֵ����
const int tnp2 = 9; // CO2�ο��¶�ȡֵ����
const int dng = 1; // ��˹���ֵ�ȡֵ���, both


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

// ����
struct Individual {

	vec1i vars; // ����
	int error; // �������
	int eqts; // ���ֵ�����
	double crowding_distance; // ӵ����

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

	int pop_size; // ��Ⱥ��ģ
	int max_gen;  // ����������
	double cross_pr; // �������
	double mutation_pr; // �������
	int mutation_num; // ��������

	int Gauss_points_num1; // H2O��˹���ֵ����
	int Gauss_points_num2; // CO2��˹���ֵ����
	int Gauss_vals_num1; // H2O��˹���ֵ�ȡ������
	int Gauss_vals_num2; // CO2��˹���ֵ�ȡ������
	int Tp_num1;  // H2O�Ĳο��¶ȸ���
	int Tp_num2;  // CO2�Ĳο��¶ȸ���
	int total_pos_num; // �ܱ����ж���λ

	int strategy_idx1; // H2O�ķ������
	int strategy_idx2; // H2O�ķ������

	vec1i val_range_map; // �洢��Ӧ����λ�õı������ͣ�1��H2O��˹���ֵ�����2��CO2��˹���ֵ�����3��H2O�ο��¶ȣ�4��CO2�ο��¶�

	int father, son; // ָʾ��ǰ��Ⱥ��������0����1

	vector<vector<Individual> > individuals; // ��Ⱥ���и���

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

		this->total_pos_num = 2 * (Gauss_points_num1 + Gauss_points_num2); // �ܱ����ж���λ��H2O��˹���ֵ���+H2O�ο��¶�+CO2��˹���ֵ���+CO2�ο��¶�+H2O���鷽��

		// ��ʼ��val_range_map
		this->val_range_map = vec1i(total_pos_num, 0);
		// H2O
		for (int i = 0; i < Gauss_points_num1; ++i) {
			val_range_map[2 * i] = 1; // H2O��˹���ֵ���
			val_range_map[2 * i + 1] = 3; // H2O�ο��¶�
		}
		// CO2
		for (int i = 0; i < Gauss_points_num2; ++i) {
			val_range_map[2 * Gauss_points_num1 + 2 * i] = 2; // CO2��˹���ֵ���
			val_range_map[2 * Gauss_points_num1 + 2 * i + 1] = 4; // CO2�ο��¶�
		}

		this->father = 0;
		this->son = 1;

		this->individuals = vector<vector<Individual> >(2, vector<Individual>(pop_size, Individual(total_pos_num)));
	}
};

// ���ݲ�ͬ�����λ���ɺ��ʵ������
inline int rand_val(VGA& data, int index) {

	const vec1i& val_range_map = data.val_range_map;

	// �������������ֵ����
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr_int_gauss1(0, data.Gauss_vals_num1 - 1); // H2O��˹���ֵ���������������
	uniform_int_distribution<int> distr_int_gauss2(0, data.Gauss_vals_num2 - 1); // CO2��˹���ֵ���������������
	uniform_int_distribution<int> distr_int_Tp1(0, data.Tp_num1 - 1); // Tp1���������������
	uniform_int_distribution<int> distr_int_Tp2(0, data.Tp_num2 - 1); // Tp2���������������

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

	// ��һ������Ⱥ���������ֵ
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
	return a.crowding_distance > b.crowding_distance; // ע���Ǵ���
}

// ӵ���������õĽṹ��
struct Front_Points {

	int index; // ��data�е�����ֵ
	int error; // �������
	int eqts; // ���ֵ�����
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
	return a.crowding_distance > b.crowding_distance; // ע���Ǵ���
}

vec1i picked_by_crowding_distance(VGA& data, vec1i& front, int num) {

	vec1i picked_nums(num);

	int& father = data.father;
	int& son = data.son;
	double delta;

	vector<Front_Points> points{}; // ��ǰfront�����и���
	for (int i = 0; i < front.size(); ++i) {
		//points.emplace_back(front[i], data.individuals[father][front[i]].error, data.individuals[father][front[i]].eqts, 0);
		points.emplace_back(data.individuals[father][front[i]], front[i]);
	}


	// ��error��С��������
	sort(points.begin(), points.end(), cmp_by_error1);
	delta = points[points.size() - 1].error - points[0].error; // ���error����Сerror�Ĳ�
	if (delta != 0) { // ��ĸ��Ϊ0

		// ��Ե����ӵ�������޴�
		points[0].crowding_distance += 1e20;
		points[points.size() - 1].crowding_distance += 1e20;

		for (int i = 1; i < points.size() - 1; ++i)
			points[i].crowding_distance += (points[i + 1].error - points[i - 1].error) / delta;
	}

	// ��eqts��С��������
	sort(points.begin(), points.end(), cmp_by_eqts1);
	delta = points[points.size() - 1].eqts - points[0].eqts; // ���eqts����Сeqts�Ĳ�
	if (delta != 0) { // ��ĸ��Ϊ0

		// ��Ե����ӵ�������޴�
		points[0].crowding_distance += 1e20;
		points[points.size() - 1].crowding_distance += 1e20;

		for (int i = 1; i < points.size() - 1; ++i)
			points[i].crowding_distance += (points[i + 1].eqts - points[i - 1].eqts) / delta;
	}

	// �������ӵ���ȴӴ�С����
	sort(points.begin(), points.end(), cmp_by_crowding_distance1);


	// ����ʣ����������ӵ�������ĸ���
	for (int i = 0; i < num; ++i) {
		picked_nums[i] = points[i].index;
	}

	return picked_nums;
}

// p֧��qΪ��
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

// ������Ӧ��
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
			// ��������Ŀ�꺯��ֵ
			individual.error = b4;  // �ø�������ֵ
			individual.eqts = 0;
			for (int ii = 0; ii < data.Gauss_points_num1; ++ii) // ע��ÿ������2
				individual.eqts += gauss_map1[individual.vars[2 * ii]];
			for (int ii = 0; ii < data.Gauss_points_num2; ++ii) // ע��ÿ������2
				individual.eqts += gauss_map2[individual.vars[2 * data.Gauss_points_num1 + 2 * ii]];
		}
	}

	return 0;
}

// ������Ҫ����ĸ�����±꼯��vec_index��ֻ������Щ�������Ӧ��
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

				b4 += (b0 + b2 + b3); // ����û�г�ww��
			}
			// ��������Ŀ�꺯��ֵ
			individual.error = b4;  // �ø�������ֵ
			individual.eqts = 0;
			for (int ii = 0; ii < data.Gauss_points_num1; ++ii) // ע��ÿ������2
			{
				individual.eqts += gauss_map1[individual.vars[2 * ii]];

			}
			for (int ii = 0; ii < data.Gauss_points_num2; ++ii) // ע��ÿ������2
			{
				individual.eqts += gauss_map2[individual.vars[2 * data.Gauss_points_num1 + 2 * ii]];
			}
		}
	}

	return 0;
}

// generate_offspring_NSGA2ר��
int non_dominated_sort_2(VGA& data, int chosen_num, int startN[][2], const vec1i& gauss_map1, const vec1i& gauss_map2, int max_mutation_repeat_times, int repeat_times) {

	// �������������ֵ����
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr_int_mutation(0, data.total_pos_num - 1); // ����ѡ��
	uniform_int_distribution<int> distr_int_mutation_gauss(0, data.Gauss_points_num1 + data.Gauss_points_num2 - 1); // ����ѡ��
	uniform_real_distribution<double> distr_double(0, 1); // 0~1���ʣ����������


	int& father = data.father;
	//int& H2O_index = data.H2O_index;
	vector<Individual>& individuals = data.individuals[father];

	sort(individuals.begin(), individuals.end(), cmp_by_error);

	//if (repeat_times > 5) {
	int err_now = individuals[0].error;
	int eqts_now = individuals[0].eqts;

	int ind_mutation = -1; // �����ڸ�����ĸ�����λ
	int mutation_bit_num = -1; // �����λ����
	vec1i vec_index{}; // ������Ҫ���¼���Ŀ�꺯���ĸ��������
	for (int i = 1; i < individuals.size(); ++i) {
		if (err_now != individuals[i].error || eqts_now != individuals[i].eqts) {
			err_now = individuals[i].error;
			eqts_now = individuals[i].eqts;
		}
		else { // ���ȡ��λ������
			//if (distr_double(eng) < 0.5) {
				//mutation_bit_num = distr_double(eng) < 0.5 ? 1 : 2; // �Ծ��ȸ��ʱ���1��2����
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


	vec2i fronts(data.pop_size, vec1i{}); // ������ǰ��
	vec1i rank(data.pop_size, 0); // �ȼ�
	vec2i S(data.pop_size, vec1i{}); // p��֧�伯��
	vec1i N(data.pop_size, 0); // p�ı�֧���

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < data.pop_size; ++i) {
		for (int j = 0; j < data.pop_size; ++j) {
			if (is_dominating(data, i, j)) { // ���i֧��j
				S[i].push_back(j); // j����i��֧�伯��
			}
			else if (is_dominating(data, j, i)) { //���i��j֧��
				N[i] += 1; //i�ı�֧�����1
			}
		}
		if (N[i] == 0)
			rank[i] = 0;
	}

	for (int i = 0; i < data.pop_size; ++i) {
		if (N[i] == 0)
			fronts[0].push_back(i);
	}

	int total_num = fronts[0].size(); // ��ѡ�е�ǰ�ص������
	int k = 0;
	vec1i Q{};
	while (fronts[k].size() != 0 && total_num < chosen_num) {
		Q.clear();
		for (int i = 0; i < fronts[k].size(); ++i) { // ��ǰ��ǰ�ص�ÿ����
			int& idx = fronts[k][i]; // ��ǰǰ�ص��ʵ�ʱ��
			//cout << "------" << idx << endl;
			for (int j = 0; j < S[idx].size(); ++j) { // ÿ����֧��ĵ�ļ���
				int& ind = S[idx][j]; // ��֧��ĵ��ʵ��������
				//cout << ind << "/ " << endl;
				N[ind] -= 1; // ��ǰǰ�ص�ȥ����
				if (N[ind] == 0) { // ���ñ�֧����ܷ��Ϊ�µ�ǰ�ص㣿
					rank[ind] = k + 1; // rank + 1
					Q.push_back(ind); // �ñ�֧��㱻������һ��ǰ�ص�
				}
			}
		}
		++k;
		fronts[k] = Q;
		total_num += Q.size(); // ��¼Ŀǰ�Ѿ�ѡ�е�ǰ�ص������
	}

	int& son = data.son;

	// ��ѡ��ӵ���������µĸ���
	int cnt = 0; // Ŀǰ��һ��������
	for (int i = 0; i < k + 1; ++i) {
		if (cnt == chosen_num)
			break;

		// Ŀǰ��һ��������������ǰ�ظ������Բ�����pop_size
		if (cnt + fronts[i].size() <= chosen_num) {
			for (int j = 0; j < fronts[i].size(); ++j) {
				int& ind = fronts[i][j];
				data.individuals[son][cnt] = data.individuals[father][ind];
				++cnt;
			}
		}
		// Ŀǰ��һ��������������ǰ�ظ���������pop_size������ӵ������ѡ
		else {
			int left_num = chosen_num - cnt; // ʣ������
			vec1i picked_nums = picked_by_crowding_distance(data, fronts[i], left_num);
			for (int j = 0; j < left_num; ++j) {
				int& ind = picked_nums[j];
				data.individuals[son][cnt] = data.individuals[father][ind];
				++cnt;
			}
		}
	}


	// ���cnt�ǲ��ǵ���pop_size
	if (cnt != chosen_num) {
		int left_num = chosen_num - cnt;
		cout << "cnt == " << cnt << ", chosen_num = " << chosen_num << endl;
		exit(-1);
	}

	return 0;
}

int generate_offspring_NSGA2(VGA& data, int GA_sgn, int beg_ind, int end_ind, int startN[][2], const vec1i& gauss_map1, const vec1i& gauss_map2, int max_mutation_repeat_times, int repeat_times) {
	double& Pc = data.cross_pr; // �������
	double& Pm = data.mutation_pr; // �������
	//int& H2O_index = data.H2O_index;
	const int& total_pos_num = data.total_pos_num;

	random_device rd;
	default_random_engine eng(rd());
	uniform_real_distribution<double> distr_double(0, 1); // 0~1���ʣ����������
	uniform_int_distribution<int> distr_int_pop(0, data.pop_size - 1); // ��Ⱥ�����������
	uniform_int_distribution<int> distr_int_cross(0, total_pos_num); // ����ѡ��
	uniform_int_distribution<int> distr_int_mutation_num(1, data.mutation_num); // ����ѡ��-total
	uniform_int_distribution<int> distr_int_mutation(0, total_pos_num - 1); // ����ѡ��-total


	int& father = data.father; // ��ǰ����
	int& son = data.son;

	int ind1, ind2; // ������ѡ��ĸ������
	double tmp1, tmp2;
	int ind_cross1, ind_cross2, ind_cross3, ind_mutation;
	double roll_val = -1;
	int mutation_point_num = -1; // ȷ��������м���

	int half_size = data.pop_size / 2;

	calc_error(data, half_size, data.pop_size, startN, gauss_map1, gauss_map2);

	// ��֧������ȡǰhalf_size�����Ӵ�
	non_dominated_sort_2(data, half_size, startN, gauss_map1, gauss_map2, max_mutation_repeat_times, repeat_times);

	// ���Ӵ���half_size�����彻�����������һ���Ӵ�
		// Pr_table�����̶��õ�����, pop_size+1��, 0��1
	double err_max = -1e10;
	vec1d Pr_table(half_size + 1, 0);
	for (int i = 0; i < half_size; ++i) {
		err_max = data.individuals[son][i].error > err_max ? data.individuals[son][i].error : err_max;
		Pr_table[i + 1] = data.individuals[son][i].error + Pr_table[i]; // 0815���£�֮ǰ��-data.errors[i] + Pr_table[i]��ɾȥ�˸���
	}

	// ������и��嶼һ���ˣ��Ѿ�����������������
	if (half_size * err_max - *(Pr_table.end() - 1) == 0) {
		return 1;
	}

	for (int i = 0; i < half_size; ++i) {
		Pr_table[i] = (i * err_max - Pr_table[i]) / (half_size * err_max - *(Pr_table.end() - 1));
	}

	// ������С��eqts����Ȩһ�����̶�
	double eqts_max = -1;
	vec1d Eq_table(half_size + 1, 0);
	for (int i = 0; i < half_size; ++i) {
		err_max = data.individuals[son][i].eqts > err_max ? data.individuals[son][i].eqts : err_max;
		Eq_table[i + 1] = data.individuals[son][i].eqts + Eq_table[i]; // 0815���£�֮ǰ��-data.errors[i] + Pr_table[i]��ɾȥ�˸���
	}
	for (int i = 0; i < half_size; ++i) {
		Eq_table[i] = (i * eqts_max - Eq_table[i]) / (half_size * eqts_max - *(Eq_table.end() - 1));
		Pr_table[i] = (Pr_table[i] + Eq_table[i]) / 2;
	}

	//vec1i break_points(3);

	// ����ͱ���
	for (int i = half_size; i < data.pop_size; ++i) { // ��elite_num֮�󣬲����µ��Ӵ���Ⱥ

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


		if (distr_double(eng) < Pc) { // ���㽻��
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
		else { // ֱ�Ӹ���
			for (int x = 0; x < total_pos_num; ++x)
				data.individuals[son][i].vars[x] = data.individuals[son][ind1].vars[x];
		}

		// �����м��������
		mutation_point_num = distr_int_mutation_num(eng);
		// ȫ�����
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

	vec1i vars; // ���б���������H2O
	double error; // �������
	int eqts; // ���ֵ�����
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

// ��ѡ�������Ÿ���
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

// �������Ӵ�
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

	vec2i fronts(data.pop_size, vec1i{}); // ������ǰ��
	vec1i rank(data.pop_size, 0); // �ȼ�
	vec2i S(data.pop_size, vec1i{}); // p��֧�伯��
	vec1i N(data.pop_size, 0); // p�ı�֧���

	int& chosen_num = solution_num;

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < data.pop_size; ++i) {
		for (int j = 0; j < data.pop_size; ++j) {
			if (is_dominating(data, i, j)) { // ���i֧��j
				S[i].push_back(j); // j����i��֧�伯��
			}
			else if (is_dominating(data, j, i)) { //���i��j֧��
				N[i] += 1; //i�ı�֧�����1
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
		for (int i = 0; i < fronts[k].size(); ++i) { // ��ǰ��ǰ�ص�ÿ����
			int& idx = fronts[k][i]; // ��ǰǰ�ص��ʵ�ʱ��
			for (int j = 0; j < S[idx].size(); ++j) { // ÿ����֧��ĵ�ļ���
				int& ind = S[idx][j]; // ��֧��ĵ��ʵ��������
				N[ind] -= 1; // ��ǰǰ�ص�ȥ����
				if (N[ind] == 0) { // ���ñ�֧����ܷ��Ϊ�µ�ǰ�ص㣿
					rank[ind] = k + 1; // rank + 1
					Q.push_back(ind); // �ñ�֧��㱻������һ��ǰ�ص�
				}
			}
		}
		++k;
		fronts[k] = Q;
	}

	vector<vector<Goal> > res{};

	// ��ѡ��ӵ���������µĸ���
	int cnt = 0; // Ŀǰ��һ��������
	vec1d errors_tmp;
	vec1i eqts_tmp;
	for (int i = 0; i < k + 1; ++i) {

		if (cnt == chosen_num)
			break;

		vector<Goal> individuals_tmp{};

		remove_repeat_point(data, fronts[i]);

		// Ŀǰ��һ��������������ǰ�ظ������Բ�����pop_size
		if (cnt + fronts[i].size() <= chosen_num) {
			for (int j = 0; j < fronts[i].size(); ++j) {
				int& ind = fronts[i][j];
				individuals_tmp.push_back(Goal(data, data.individuals[data.father][ind], gauss_map1, gauss_map2));
				++cnt;
			}
		}
		// Ŀǰ��һ��������������ǰ�ظ���������pop_size������ӵ������ѡ
		else {
			int left_num = chosen_num - cnt; // ʣ������
			vec1i picked_nums = picked_by_crowding_distance(data, fronts[i], left_num);
			for (int j = 0; j < left_num; ++j) {
				int& ind = picked_nums[j];
				individuals_tmp.push_back(Goal(data, data.individuals[data.father][ind], gauss_map1, gauss_map2));
				++cnt;
			}
		}

		res.push_back(individuals_tmp);

	}

	// ���cnt�ǲ��ǵ���pop_size
	//if (cnt != chosen_num) {
	//	cout << "[ERROR] cnt == " << cnt << ", chosen_num = " << chosen_num << endl;
	//	exit(-1);
	//}

	return res;
}


int NSGA2_one(int strategy_idx1, int strategy_idx2, int pop_size, int max_gen, string method_name, int startN[][2], Goal& err_now) {
	// �Ŵ��㷨begin

	double total_it = 0; // ��¼�ܵ�������
	double total_time = 0; // ��¼��ʱ��

	// ��ʱ��ʼ
	auto start = chrono::steady_clock::now();

	double error_min = 1e20;

	const int vars_num = GN1 + GN2; // ����������ˮ����GN1λ+������̼GN2λ+ˮ�������鷽��1λ
	const int max_muation_bit = int(vars_num / 2); // ������λ��
	constexpr double cross_pr = 0.8; // ���㽻�����


	// ����������
	VGA data(pop_size, max_gen, GN1, GN2, tng1, tng2, strategy_idx1, strategy_idx2, tnp1, tnp2);

	int elite_num = int(0.05 * pop_size); // ��Ӣ��������

	// ��ʼ����Ⱥ
	init_GA(data);

	vec1i gauss_map1{}, gauss_map2{};
	for (int i = 0; i < tng1; ++i)
		gauss_map1.push_back(ngls1 + dng * i);
	for (int i = 0; i < tng2; ++i)
		gauss_map2.push_back(ngls2 + dng * i);

	calc_error(data, 0, pop_size, startN, gauss_map1, gauss_map2);

	Goal err_best{}; // �ճ�ʼ�������ڲ���ʹ��

	int repeat_times = 0;

	char folder_name[256];
	snprintf(folder_name, sizeof(folder_name) - 1, ".\\%s", method_name);
	mkdir(folder_name);

	char folder_name_detail[256];
	snprintf(folder_name_detail, sizeof(folder_name_detail) - 1, "%s\\detail-%d", folder_name, pop_size);
	mkdir(folder_name_detail);

	FILE* fp;
	int removed_it = 50; // ȥ���ظ�����ѭ������
	int GA_sgn = 3; // 1��GA�� 2��NSGA2

	int max_repeat_times = 50; // �ﵽ���Ž�ͣ�͵������������ѭ��
	int max_mutation_repeat_times = int(0.2 * max_repeat_times);
	int max_mutation_repeat_times_2 = int(0.1 * max_repeat_times); // ������������Ž�ͣ�͵����������������ķ�ĸ
	int max_mutation_repeat_times_1 = int(0.5 * max_repeat_times); // ������������Ž�ͣ�͵����������������ķ�ĸ
	int change_method_times = int(0.1 * max_repeat_times); // ÿ�������λ�GA����
	int only_NSGA2_times = int(0.6 * max_repeat_times); // �����õ���������ֻʹ��NSGA2

	// �Ŵ��㷨����begin
	for (int it = 0; it < max_gen; ++it, ++repeat_times) {

		printf("%s generation %d, %d, %.3f...\r", method_name, it, repeat_times, data.mutation_pr);


		// �����Ӵ���������Ӧ�ȡ�ѡ�񡢽��桢����
		generate_offspring_NSGA2(data, GA_sgn, 0, pop_size, startN, gauss_map1, gauss_map2, max_mutation_repeat_times, repeat_times);

		// ��ǰ�����Ž�
		err_best = find_goal(data, gauss_map1, gauss_map2);

		// ����и��Ž�
		if (err_best.error < error_min) {

			// ����̨���
			printf("\n%d || %d || %.2f %d || %d %d || ", strategy_idx2, it, err_best.error, err_best.eqts, err_best.H2O_group_num, strategy_idx2);
			for (int i = 0; i < err_best.total_pos_num; ++i) {
				printf("%d ", err_best.vars[i]);
			}
			printf("\n");
			error_min = err_best.error;

			// �ļ����
			char name[256];
			snprintf(name, sizeof(name) - 1, "%s\\res-%d-%d.txt", folder_name_detail, strategy_idx1, strategy_idx2);
			fp = fopen(name, "a");
			fprintf(fp, "%d %.4f %d %d %d || ", it, err_best.error, err_best.eqts, err_best.H2O_group_num, strategy_idx2);
			for (int i = 0; i < err_best.total_pos_num; ++i) {
				fprintf(fp, "%d ", err_best.vars[i]);
			}
			fprintf(fp, "\n");
			fclose(fp);

			// ���Ž�ͣ�͵Ĵ�������
			repeat_times = 0;
			//repeat_times = int(0.1 * max_mutation_repeat_times);
		}

		// ������Ž�ͣ�͵Ĵ�������500��������ѭ���������һ�Σ�д�ļ�����������
		if (repeat_times > max_repeat_times || it == max_gen - 1) {

			// �ļ����
			char name[256];
			snprintf(name, sizeof(name) - 1, "%s\\res-%d-%d.txt", folder_name_detail, strategy_idx1, strategy_idx2);
			fp = fopen(name, "a");
			fprintf(fp, "%d %.4f %d %d %d || ", it, err_best.error, err_best.eqts, err_best.H2O_group_num, strategy_idx2);
			for (int i = 0; i < err_best.total_pos_num; ++i) {
				fprintf(fp, "%d ", err_best.vars[i]);
			}
			fprintf(fp, "\n");
			fclose(fp);

			total_it = it; // ��¼�ܵ�������

			break; // ��������
		}


		// ������Ž�ͣ�͵Ĵ�������100�������ӱ����λ����������ǿ��Ķ����ԣ�����������Χ
		else if (repeat_times > 20) {
			data.mutation_num = 2;
			data.mutation_pr = 0.2;
		}
		// �������Ž�ͣ�͵Ĵ��������࣬�������������1
		else {
			data.mutation_num = 1;
			data.mutation_pr = 0.05;
		}
		// ��ת�����Ӵ���ʶ
		swap_father_son_sgn(data);

	}

	// �Ŵ��㷨����end
	printf("\n");

	// д��ǰ������Ž�
	{

		// д��ǰ���Ž�2
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

		// дparetoǰ�ؽ����ÿ��zone����һ��ǰ��
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

		// дparetoǰ����ϸ�����ÿ��zone����һ��ǰ��
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

	// ��ʱ����;
	auto end = chrono::steady_clock::now();
	auto time_diff = end - start;
	auto duration = chrono::duration_cast<chrono::seconds>(time_diff);
	cout << "Operation cost : " << duration.count() << "s" << endl;

	total_time = chrono::duration_cast<chrono::milliseconds>(time_diff).count() / 1e3; // ��¼��ʱ��(s)

	err_now = err_best;

	return 0;
	//return vec1d{ total_it, total_time, total_time / total_it };
}

