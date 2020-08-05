/*
	Calculates the effective reproduction number and outbreak control
	with using 3 interventions:
		1) case isolation
		2) contact tracing
		3) mask wearing

	[ based on the work of Hellewell J., Abbott S., Gimma A., et al. https://github.com/cmmid/ringbp ]

*/

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <limits>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

using namespace std;

double Inf = numeric_limits<double>::infinity();


void delayDistribution(double x_start, double x_end, double step) // Weibull distribution: short delay = {shape=1.651524; scale=4.287786}; long delay = {shape=2.305172; scale=9.483875}
{
	// Parameters for short delay (for the late stages of the SARS 2003 outbreak in Singapore)
	double shape_short = 1.651524;
	double scale_short = 4.287786;
	boost::math::weibull_distribution<double> weibull_short(shape_short, scale_short);

	// Parameters for long delay (for early phase of the COVID-19 outbreak in Wuhan)
	double shape_long = 2.305172;
	double scale_long = 9.483875;
	boost::math::weibull_distribution<double> weibull_long(shape_long, scale_long);

	ofstream fout;
	fout.open("delayDist.csv");
	fout << "Time since infection (days)" << "," << "PDF short delay" << "," << "PDF long delay" << "\n";
	double x = x_start;
	double density_short = 0.0;
	double density_long = 0.0;


	while (x < x_end)
	{
		density_short = pdf(weibull_short, x);
		density_long = pdf(weibull_long, x);

		fout << setprecision(8) << x << "," << density_short << "," << density_long << "\n";
		x = x + step;
	}

	fout.close();
}


void incubationDistribution(double x_start, double x_end, double step) // Weibull distribution: shape = 2.322737; scale = 6.492272
{
	double shape = 2.322737;
	double scale = 6.492272;
	boost::math::weibull_distribution<double> weibull(shape, scale);

	cout << mean(weibull);
	cin.get();

	ofstream fout;
	fout.open("incubationDist.csv");
	fout << "Time since infection (days)" << "," << "PDF" << "\n";
	double x = x_start;
	double density = 0.0;

	while (x < x_end)
	{
		density = pdf(weibull, x);

		fout << setprecision(8) << x << "," << density << "\n";
		x = x + step;
	}

	fout.close();
}


void serialIntervalDistribution(double x_start, double x_end, double step) // Skew normal ditribution: location = incubation period (eg, 5); scale = 2; shape = 30 | 1.95 | 
{
	double location = 5.0;		// mean if normal --> the incubation period of the case (5.0 - just example)
	double scale = 2.0;			// width; standart deviation if normal
	double shape_1p = 30.0;		// skew; the distribution is right skewed if shape > 0 and is left skewed if shape < 0; shape = 0 => normal distribution --> proportion of transmission <1%
	double shape_15p = 1.95;	// skew --> proportion of transmission before symptoms 15%
	double shape_30p = 0.7;		// skew --> proportion of transmission before symptoms 30%
	double shape_40p = 0.324;	// skew --> proportion of transmission before symptoms 40%

	double location_asym = 19; // mean = meadian for asymptomatic cases
	double std = 11.0 / 1.349; // IQR = [15; 26];	1.349 - value from IQR for Normal distribution
	x_end = 30;
	boost::math::skew_normal_distribution<double> normal_asym(location_asym, std, 0.0);

	boost::math::skew_normal_distribution<double> skew_normal_1p(location, scale, shape_1p);
	boost::math::skew_normal_distribution<double> skew_normal_15p(location, scale, shape_15p);
	boost::math::skew_normal_distribution<double> skew_normal_30p(location, scale, shape_30p);
	boost::math::skew_normal_distribution<double> skew_normal_40p(location, scale, shape_40p);
	boost::math::skew_normal_distribution<double> normal(location, scale, 0.0);
	// play with Gamma distibution
	double disp_param = 0.16;
	double mean_value = 3.5; // initial R0

	double prob_success1 = disp_param / (disp_param + mean_value);
	double prob_success2 = 0.13 / (0.13 + mean_value);

	double gamma_scale1 = (1.0 - prob_success1) / prob_success1;
	double gamma_scale2 = (1.0 - prob_success2) / prob_success2;
	if (gamma_scale1 == 0.0 || gamma_scale2 == 0.0)
	{
		gamma_scale1 = 0.0000000000001;
		gamma_scale2 = 0.0000000000001;
	}
	//gamma_distribution<> gamma(disp_param, gamma_scale);
	boost::math::gamma_distribution<> gamma1(disp_param, gamma_scale1);
	boost::math::gamma_distribution<> gamma2(0.13, gamma_scale2);

	boost::math::poisson_distribution<> poisson1(mean(gamma1));
	boost::math::poisson_distribution<> poisson2(mean(gamma2));

	cout << "standart case:" << endl;
	cout << "Mean Gamma(shape=" << disp_param << "; scale=" << gamma_scale1 << ") = " << mean(gamma1) << "; probability of success = " << prob_success1 << endl;
	cout << "Mean Poisson(rate=" << mean(gamma1) << ") = " << mean(poisson1) << endl;
	cout << "aka add a wearing of masks:" << endl;
	cout << "Mean Gamma(shape=" << disp_param << "; scale=" << gamma_scale2 << ") = " << mean(gamma2) << "; probability of success = " << prob_success2 << endl;
	cout << "Mean Poisson(rate=" << mean(gamma2) << ") = " << mean(poisson2) << endl;
	cin.get();
	//poisson_distribution<> poisson(gamma(rng));


	//case_data.new_cases[i] = poisson(rng);

	// check 
	//cout << "The probability that random X will less than or equal to incubation period: " << cdf(normal, location) << endl;
	//cin.get();

	ofstream fout;
	fout.open("serialIntervalDist.csv");
	//fout.open("negBinomDist.csv");
	//fout << "Mean value" << "," << setprecision(2) << location << "\n";
	//fout << "Time since infection (days)" << "," << "PDF <1%" << "," << "PDF 15%" << "," << "PDF 30%" << "," << "PDF 40%" << "," << "PDF 50% (normal)" << "\n";
	fout << "Time since infection (days)" << "," << "PDF 40%" << "," << "PDF normal asym" << "\n";
	//fout << "Time since infection (days)" << "," << "PDF standart" << "," << "PDF with masks" << "\n";
	double x = x_start;
	double density_1p = 0.0;
	double density_15p = 0.0;
	double density_30p = 0.0;
	double density_40p = 0.0;
	double density1 = 0.0;
	double density2 = 0.0;
	double density = 0.0;

	while (x < x_end)
	{
		density_1p = pdf(skew_normal_1p, x);
		density_15p = pdf(skew_normal_15p, x);
		density_30p = pdf(skew_normal_30p, x);
		density_40p = pdf(skew_normal_40p, x);
		density = pdf(normal_asym, x);
		density1 = pdf(poisson1, x);
		density2 = pdf(poisson2, x);

		//fout << setprecision(8) << x << "," << density_1p << "," << density_15p << "," << density_30p << "," << density_40p << "," << density << "\n";
		fout << setprecision(8) << x << "," << density_40p << "," << density << "\n";
		//fout << setprecision(8) << x << "," << density1 << "," << density2 << "\n";
		x = x + step;
	}

	fout.close();
}


struct Case_data
{
	int array_size;
	vector<bool> mask; // wearing mask or not
	vector<double> exposure; //time of getting infection
	vector<bool> asym; // asymptomatic or not
	vector<int> case_id;
	vector<int> infector;
	vector<bool> missed; // traced or not
	vector<double> onset; //time of symptoms onset;
	vector<int> new_cases;
	vector<double> isolated_time;
	vector<bool> isolated;
};


/* Structure for output from outbreak_step() function */
struct Outbreak_step_output
{
	struct Case_data case_data;
	double effective_R0;
	int cases_in_gen;				// number of new cases
	int cases_symptomatic;
	int cases_asymptomatic;
};


/* Structure for output from outbreak_model() function (effective R0 and cases per generation*/
struct Model_output
{
	vector<double> effective_R0;
	double mean_effective_R0;
	int total_cases;
	int total_symptomatic;
	int total_asymptomatic;
	bool controlled;
	vector<int> cases_in_gen;
	vector<double> latest_onset;
	vector<int> weekly_cases;
};

// for checking
void show_data(struct Case_data case_data)
{
	cout << "exposure | asym | case_id | infector | missed | onset | new_cases | isolated_time | isolated " << endl;
	for (int i = 0; i < case_data.array_size; i++)
	{
		cout << setw(9) << setprecision(3) << case_data.exposure[i] << "|";
		cout << setw(6) << boolalpha << case_data.asym[i] << "|";			// boolalpha --> print "true/false" instead of "1/0"
		cout << setw(9) << case_data.case_id[i] << "|";
		cout << setw(10) << case_data.infector[i] << "|";
		cout << setw(8) << case_data.missed[i] << "|";
		cout << setw(7) << setprecision(3) << case_data.onset[i] << "|";
		cout << setw(11) << case_data.new_cases[i] << "|";
		cout << setw(15) << setprecision(3) << case_data.isolated_time[i] << "|";
		cout << setw(9) << case_data.isolated[i] << " \n";
	}

}

// for checking
void show_data(struct Outbreak_step_output out)
{
	cout << "exposure | asym | case_id | infector | missed | onset | new_cases | isolated_time | isolated " << endl;
	for (int i = 0; i < out.case_data.array_size; i++)
	{
		cout << setw(9) << setprecision(3) << out.case_data.exposure[i] << "|";
		cout << setw(6) << out.case_data.asym[i] << "|";
		cout << setw(9) << out.case_data.case_id[i] << "|";
		cout << setw(10) << out.case_data.infector[i] << "|";
		cout << setw(8) << out.case_data.missed[i] << "|";
		cout << setw(7) << setprecision(3) << out.case_data.onset[i] << "|";
		cout << setw(11) << out.case_data.new_cases[i] << "|";
		cout << setw(15) << setprecision(3) << out.case_data.isolated_time[i] << "|";
		cout << setw(9) << out.case_data.isolated[i] << " \n";
	}
	cout << setprecision(5) << "Effective R0 = " << out.effective_R0 << endl;
	cout << "Cases in generation = " << out.cases_in_gen << endl;
}


vector<double> get_inc_period(int n)
{
	vector<double> inc_period(n);
	double shape = 2.322737;
	double scale = 6.492272;
	random_device dev;
	default_random_engine rng(dev()); // random numbers generator
	weibull_distribution<double> weibull(shape, scale);

	for (int i = 0; i < n; i++)
	{
		inc_period[i] = weibull(rng);
	}

	return inc_period;
}


double* get_delay(int n, double shape, double scale)
{
	double* delay = new double[n];
	random_device dev;
	default_random_engine rng(dev()); // random numbers generator
	weibull_distribution<double> weibull(shape, scale);

	for (int i = 0; i < n; i++)
	{
		delay[i] = weibull(rng);
	}

	return delay;
}


double get_serial_interval(double mean, double SD, double k) // for symptomatic cases
{
	random_device dev;
	default_random_engine rng(dev()); // random numbers generator
	normal_distribution<double> normal(0, 1); // standart normal distribution
	double u0 = normal(rng);
	double v = normal(rng);
	double u1, delta, x;

	if (k == 0) // if no skewness
		x = mean + SD * u0;
	else
	{
		// add skewness
		delta = k / sqrt(1 + k * k);
		u1 = delta * u0 + sqrt(1.0 - delta * delta) * v;
		if (u0 < 0)
			u1 = -u1;
		x = mean + SD * u1;
	}

	if (x < 1.0)
		return 1.0;
	return x;
}

double get_asym_serial_interval(double mean, double std) // for asymptomatic cases
{
	random_device dev;
	default_random_engine rng(dev()); // random numbers generator
	normal_distribution<double> normal(mean, std); // standart normal distribution

	double x = normal(rng);

	if (x < 1.0)
		return 1.0;

	return x;
}

struct Case_data outbreak_setup(double num_cases, double num_cases_asym, double delay_shape, double delay_scale, double prob_asym, double prob_mask)
{
	random_device dev;
	default_random_engine rng(dev()); // random numbers generator
	bernoulli_distribution bernoulli(prob_asym);
	bernoulli_distribution bernoulli_mask(prob_mask);
	double delay = *get_delay(1, delay_shape, delay_scale);
	struct Case_data case_data;
	int total_cases = num_cases + num_cases_asym;
	//int total_cases = num_cases;

	case_data.array_size = total_cases;
	case_data.mask = vector<bool>(total_cases);
	case_data.exposure = vector<double>(total_cases);
	case_data.asym = vector<bool>(total_cases);
	case_data.case_id = vector<int>(total_cases);
	case_data.infector = vector<int>(total_cases);
	case_data.missed = vector<bool>(total_cases);
	case_data.onset = vector<double>(total_cases);
	case_data.new_cases = vector<int>(total_cases);
	case_data.isolated_time = vector<double>(total_cases);
	case_data.isolated = vector<bool>(total_cases);

	case_data.onset = get_inc_period(total_cases);
	for (int i = 0; i < num_cases; i++)
	{
		case_data.mask[i] = bernoulli_mask(rng);
		case_data.case_id[i] = i + 1;
		case_data.exposure[i] = 0.0;
		//case_data.asym[i] = bernoulli(rng);
		case_data.asym[i] = false;
		case_data.infector[i] = 0.0;
		case_data.missed[i] = true;
		case_data.new_cases[i] = -1;
		case_data.isolated_time[i] = case_data.onset[i] + delay;
		case_data.isolated[i] = false;

		//if (case_data.asym[i])
			//case_data.isolated_time[i] = Inf;
	}

	for (int i = num_cases; i < case_data.array_size; i++)
	{
		case_data.mask[i] = bernoulli_mask(rng);
		case_data.case_id[i] = i + 1;
		case_data.exposure[i] = 0.0;
		case_data.asym[i] = true;
		case_data.infector[i] = 0.0;
		case_data.missed[i] = true;
		case_data.new_cases[i] = -1;
		case_data.isolated_time[i] = Inf;
		case_data.isolated[i] = false;
	}


	return case_data;
}



struct Outbreak_step_output outbreak_step(Case_data case_data, double R0_iso, double R0_com, double disp_iso, double disp_com, double k, double delay_shape, double delay_scale,
	double prob_ascertain, double prob_asym, double prob_asym_asym, double inf_asym, double quarantine, double prob_mask, double eff_mask, bool after_symptom, bool all)
{
	struct Outbreak_step_output out;

	double prob_success_com = 0.0;
	double prob_success_mask = 0.0;
	double disp_param = 0.0;
	double mean_com = 0.0;
	double mean_mask = 0.0;
	double gamma_scale_com = 0.0;
	double gamma_scale_mask = 0.0;
	random_device dev;
	default_random_engine rng(dev()); // random numbers generator
	bernoulli_distribution bernoulli(prob_asym); // asym from symptomatic
	bernoulli_distribution bernoulli_asym(prob_asym_asym); // asym from asym
	bernoulli_distribution bernoulli_2(1.0 - prob_ascertain);
	bernoulli_distribution bernoulli_mask_wearing(prob_mask);
	bernoulli_distribution bernoulli_mask_protect(eff_mask);

	int total_cases_mask = 0;
	vector<int> new_cases_mask(case_data.array_size); // number of new cases with using reduced R0 --> infector wears mask
	vector<bool> mask_wearing;		// who from new cases will be in mask

	/* For each case in case_data, draw new_cases from a negative binomial distribution with an R0 and dispersion dependent on if isolated = TRUE */
	for (int i = 0; i < case_data.array_size; i++)
	{
		if (case_data.isolated[i])
		{
			disp_param = disp_iso;
			mean_com = R0_iso;
			mean_mask = R0_iso;
		}
		else
		{
			disp_param = disp_com;
			mean_com = R0_com;
			mean_mask = R0_com * (1.0 - eff_mask);

			if (case_data.mask[i] && !after_symptom) // if people start to wear mask after exposure
				mean_com = R0_com * (1.0 - eff_mask);

			if (case_data.asym[i]) // if asymptomatic case
			{
				mean_com = R0_com * inf_asym;
				mean_mask = mean_com * (1.0 - eff_mask);
			}
		}

		prob_success_com = disp_param / (disp_param + mean_com);

		prob_success_mask = disp_param / (disp_param + mean_mask); // only use for wearing after symptoms 

		//negative_binomial_distribution<> neg_binomial(disp_param, prob_success);
		// negative binomial = combination of Gamma and Poisson distributions
		gamma_scale_com = (1.0 - prob_success_com) / prob_success_com;
		gamma_scale_mask = (1.0 - prob_success_mask) / prob_success_mask;
		if (gamma_scale_com == 0.0)
			gamma_scale_com = 0.0000000000001;
		if (gamma_scale_mask == 0.0)
			gamma_scale_mask = 0.0000000000001;

		gamma_distribution<> gamma_com(disp_param, gamma_scale_com);
		poisson_distribution<> poisson_com(gamma_com(rng));

		gamma_distribution<> gamma_mask(disp_param, gamma_scale_mask);
		poisson_distribution<> poisson_mask(gamma_mask(rng));

		//case_data.new_cases[i] = neg_binomial(rng);
		// 1 -- generate new cases with initial R0 (community)
		case_data.new_cases[i] = poisson_com(rng);
		int cases = case_data.new_cases[i];
		// if wearing all people -- infected and suspected
		if (all)
		{
			for (int j = 0; j < case_data.new_cases[i]; j++)
			{
				bool in_mask = bernoulli_mask_wearing(rng); // define - this person wears mask or not
				bool luck = bernoulli_mask_protect(rng);	// define - would be this person protect by his mask
				if (in_mask && luck)
					cases = cases - 1;
				else
					mask_wearing.push_back(in_mask);
			}
			case_data.new_cases[i] = cases;
		}

		// 2 -- generate new cases with reduced R0 (mask)
		new_cases_mask[i] = poisson_mask(rng);
	}



	vector<double> exposure_time_all; // get_serial_interval(new_case_data.onset[i], 2.0, k);
	if (after_symptom) // if people start to wear mask after symptom (not all time after exposure)
	{
		/* Predict exposure time for new cases */
		int precit_new_cases_com = 0;
		int precit_new_cases_mask = 0;
		int presymptom_com = 0;
		int aftersymptom_com = 0;
		int aftersymptom_mask = 0;

		vector<int> total_mixed_cases(case_data.array_size);

		vector<double> exposure_time_after_com; // get_serial_interval(new_case_data.onset[i], 2.0, k);
		vector<double> exposure_time_after_mask; // get_serial_interval(new_case_data.onset[i], 2.0, k);

		double T = 0.0;
		for (int i = 0; i < case_data.array_size; i++)
		{
			precit_new_cases_com = case_data.new_cases[i];
			precit_new_cases_mask = new_cases_mask[i];

			// # presymptomatic cases without masks  
			presymptom_com = 0;
			exposure_time_after_com.clear();
			for (int j = 0; j < precit_new_cases_com; j++)
			{
				T = get_serial_interval(case_data.onset[i], 2.0, k);
				if (T < case_data.onset[i])
				{
					presymptom_com++;
					exposure_time_all.push_back(T);
				}
				else
					exposure_time_after_com.push_back(T);
			}

			// # aftersymptom cases with masks
			aftersymptom_mask = 0;
			exposure_time_after_mask.clear();
			if (case_data.mask[i])
			{
				for (int j = 0; j < precit_new_cases_mask; j++)
				{
					T = get_serial_interval(case_data.onset[i], 2.0, k);
					if (T > case_data.onset[i])
					{
						aftersymptom_mask++;
						exposure_time_after_mask.push_back(T);
					}
				}

				// # total new cases
				total_mixed_cases[i] = presymptom_com + aftersymptom_mask;
				exposure_time_all.insert(exposure_time_all.end(), exposure_time_after_mask.begin(), exposure_time_after_mask.end());
			}
			else
			{
				// # total new cases
				total_mixed_cases[i] = precit_new_cases_com;
				exposure_time_all.insert(exposure_time_all.end(), exposure_time_after_com.begin(), exposure_time_after_com.end());
			}

		}

		case_data.new_cases.clear();
		case_data.new_cases = total_mixed_cases;
	}


	/* Select cases that have generated any new cases */
	struct Case_data new_case_data;
	int total_new_cases = 0;
	int size = 0;
	for (int i = 0; i < case_data.array_size; i++)
	{
		if (case_data.new_cases[i] > 0)
		{
			total_new_cases = total_new_cases + case_data.new_cases[i];
			size++;
		}
	}
	// Allocate memory for new_case_data
	new_case_data.array_size = size;

	new_case_data.mask = vector<bool>(size);
	new_case_data.exposure = vector<double>(size);
	new_case_data.asym = vector<bool>(size);
	new_case_data.case_id = vector<int>(size);
	new_case_data.infector = vector<int>(size);
	new_case_data.missed = vector<bool>(size);
	new_case_data.onset = vector<double>(size);
	new_case_data.new_cases = vector<int>(size);
	new_case_data.isolated_time = vector<double>(size);
	new_case_data.isolated = vector<bool>(size);
	// Values setup for new_cases_data
	int j = 0;
	for (int i = 0; i < case_data.array_size; i++)
	{
		if (case_data.new_cases[i] > 0)
		{
			new_case_data.mask[j] = case_data.mask[i];
			new_case_data.exposure[j] = case_data.exposure[i];
			new_case_data.asym[j] = case_data.asym[i];
			new_case_data.case_id[j] = case_data.case_id[i];
			new_case_data.infector[j] = case_data.infector[i];
			new_case_data.isolated[j] = case_data.isolated[i];
			new_case_data.isolated_time[j] = case_data.isolated_time[i];
			new_case_data.missed[j] = case_data.missed[i];
			new_case_data.new_cases[j] = case_data.new_cases[i];
			new_case_data.onset[j] = case_data.onset[i];
			j++;
		}
	}

	/* If no new cases drawn, outbreak is over so return case_data */
	if (total_new_cases == 0)
	{
		out.cases_asymptomatic = 0;
		out.cases_symptomatic = 0;
		for (int i = 0; i < case_data.array_size; i++)
		{
			case_data.isolated[i] = true;
			if (case_data.asym[i])
				out.cases_asymptomatic++;
			else
				out.cases_symptomatic++;
		}

		out.case_data = case_data;
		out.effective_R0 = 0.0;
		out.cases_in_gen = 0;

		return out;
	}

	/* If we have new cases */
	double mean_asym_serial = 19.0;
	double std_asym_serial = 11.0 / 1.349;

	vector<double> inc_samples = get_inc_period(total_new_cases);		// sample from the incubation period for each new person

	vector<bool> mask(total_new_cases);
	vector<double> exposure(total_new_cases); //  time when new cases were exposed, a draw from serial interval based on infector's onset
	vector<int> infector(total_new_cases); // infector of each new person
	vector<double> infector_iso_time(total_new_cases); // when infector was isolated
	vector<bool> infector_asym(total_new_cases); // if infector is asymptomatic
	vector<bool> asym(total_new_cases); // draws a sample to see if this person is asymptomatic
	vector<bool> missed(total_new_cases); // draws a sample to see if this person is traced
	vector<bool> isolated(total_new_cases);
	vector<int> new_cases(total_new_cases);

	int index = 0;
	for (int i = 0; i < new_case_data.array_size; i++)
	{
		int number_cases = new_case_data.new_cases[i];

		for (int j = 0; j < number_cases; j++)
		{
			infector_asym[index] = new_case_data.asym[i];

			if (all)
				mask[index] = mask_wearing[index];
			else
				mask[index] = bernoulli_mask_wearing(rng); // define - this person wears mask or not

			if (after_symptom) // don't use this scenario for asymptomatic cases
				exposure[index] = exposure_time_all[index]; // get_serial_interval(new_case_data.onset[i], 2.0, k);
			else
			{
				if (!infector_asym[index]) // if parent is symptomatic
					exposure[index] = get_serial_interval(new_case_data.onset[i], 2.0, k);
				else
					exposure[index] = get_asym_serial_interval(mean_asym_serial, std_asym_serial);
			}

			if (!infector_asym[index]) // if parent is symptomatic
				asym[index] = bernoulli(rng);
			else
				asym[index] = bernoulli_asym(rng);

			infector[index] = new_case_data.case_id[i];
			infector_iso_time[index] = new_case_data.isolated_time[i];
			missed[index] = bernoulli_2(rng);
			if (asym[index])
				missed[index] = true;
			isolated[index] = false;
			new_cases[index] = -1;
			index++;
		}
	}

	// filter out new cases prevented by isolation (numer of these cases = new_size)
	int new_size = 0;
	for (int i = 0; i < total_new_cases; i++)
	{
		if (exposure[i] < infector_iso_time[i])
			new_size++;
	}

	vector<double> inc_samples2(new_size);
	vector<bool> mask2(new_size);
	vector<double> exposure2(new_size);
	vector<int> infector2(new_size);
	vector<double> infector_iso_time2(new_size);
	vector<bool> infector_asym2(new_size);
	vector<bool> asym2(new_size);
	vector<bool> missed2(new_size);
	vector<bool> isolated2(new_size);
	vector<int> new_cases2(new_size);
	vector<double> onset(new_size);

	index = 0;
	for (int i = 0; i < total_new_cases; i++)
	{
		if (exposure[i] < infector_iso_time[i])
		{
			mask2[index] = mask[i];
			inc_samples2[index] = inc_samples[i];
			exposure2[index] = exposure[i];
			infector2[index] = infector[i];
			infector_iso_time2[index] = infector_iso_time[i];
			infector_asym2[index] = infector_asym[i];
			asym2[index] = asym[i];
			missed2[index] = missed[i];
			isolated2[index] = isolated[i];
			new_cases2[index] = new_cases[i];
			onset[index] = exposure[i] + inc_samples[i];
			index++;
		}
	}

	// Addition operations
	// Idetify isolated time for each case
	double d = *get_delay(1, delay_shape, delay_scale);
	vector<double> isolated_time(new_size);
	for (int i = 0; i < new_size; i++)
	{
		// Cases whose parents are asymptomatic are automatically missed
		if (infector_asym2[i])
			missed2[i] = true;

		// If you are asymptomatic, your isolation time is Inf
		if (asym2[i])
			isolated_time[i] = Inf;
		else
		{
			// If you are not asymptomatic, but you are missed, you are isolated at your symptom onset
			if (missed2[i])
			{
				isolated_time[i] = onset[i] + d;
			}
			else
			{
				// If you are not asymptomatic and you are traced, your isolation time is following
				if (!quarantine)
				{
					isolated_time[i] = min(onset[i] + d, max(onset[i], infector_iso_time2[i]));
				}
				else
					isolated_time[i] = infector_iso_time2[i];
			}
		}

	}

	int not_isolated_cases = 0;
	for (int i = 0; i < case_data.array_size; i++)
	{
		if (!case_data.isolated[i])
		{
			not_isolated_cases++;
			// Everyone in case_data so far has had their chance to infect and are therefore considered isolated
			//if (!case_data.asym[i])
			case_data.isolated[i] = true;
		}
	}
	// Output data
	int N = case_data.array_size + new_size;	// total cases after this step
	out.cases_symptomatic = 0;					// total symptomatic cases after this step
	out.cases_asymptomatic = 0;					// total symptomatic cases after this step
	out.case_data.array_size = N;
	out.cases_in_gen = new_size;
	out.effective_R0 = 1.0 * new_size / not_isolated_cases;
	out.case_data.mask = vector<bool>(N);
	out.case_data.exposure = vector<double>(N);
	out.case_data.asym = vector<bool>(N);
	out.case_data.case_id = vector<int>(N);
	out.case_data.infector = vector<int>(N);
	out.case_data.missed = vector<bool>(N);
	out.case_data.onset = vector<double>(N);
	out.case_data.new_cases = vector<int>(N);
	out.case_data.isolated_time = vector<double>(N);
	out.case_data.isolated = vector<bool>(N);


	for (int i = 0; i < case_data.array_size; i++)
	{
		out.case_data.mask[i] = case_data.mask[i];
		out.case_data.exposure[i] = case_data.exposure[i];
		out.case_data.asym[i] = case_data.asym[i];
		out.case_data.case_id[i] = case_data.case_id[i];
		out.case_data.infector[i] = case_data.infector[i];
		out.case_data.missed[i] = case_data.missed[i];
		out.case_data.onset[i] = case_data.onset[i];
		out.case_data.new_cases[i] = case_data.new_cases[i];
		out.case_data.isolated_time[i] = case_data.isolated_time[i];
		out.case_data.isolated[i] = case_data.isolated[i];
		if (case_data.asym[i])
			out.cases_asymptomatic++;
	}
	index = 0;
	for (int i = case_data.array_size; i < N; i++)
	{
		out.case_data.case_id[i] = out.case_data.case_id[i - 1] + 1;
		out.case_data.mask[i] = mask2[index];
		out.case_data.exposure[i] = exposure2[index];
		out.case_data.asym[i] = asym2[index];
		out.case_data.infector[i] = infector2[index];
		out.case_data.missed[i] = missed2[index];
		out.case_data.onset[i] = onset[index];
		out.case_data.new_cases[i] = new_cases2[index];
		out.case_data.isolated_time[i] = isolated_time[index];
		out.case_data.isolated[i] = isolated2[index];
		if (asym2[index])
			out.cases_asymptomatic++;
		index++;
	}

	out.cases_symptomatic = N - out.cases_asymptomatic;

	return out;
}



struct Model_output outbreak_model(double num_initial_cases, double num_cases_asym, int max_days, double max_cases, double R0_iso, double R0_com,
	double disp_iso, double disp_com, double k, double delay_shape, double delay_scale,
	double prob_ascertain, double prob_asym, double prob_asym_asym, double inf_asym, double quarantine, double prob_mask, double eff_mask, bool after_symptom, bool all)
{
	double total_cases = num_initial_cases;
	double latest_onset = 0.0;
	bool extinct = false;

	struct Case_data case_data = outbreak_setup(num_initial_cases, num_cases_asym, delay_shape, delay_scale, prob_asym, prob_mask);
	//show_data(case_data);
	struct Outbreak_step_output out;
	struct Model_output model;
	int step = 1;
	while (latest_onset < max_days && total_cases < max_cases && !extinct)
	{
		out = outbreak_step(case_data, R0_iso, R0_com, disp_iso, disp_com, k,
			delay_shape, delay_scale, prob_ascertain, prob_asym, prob_asym_asym, inf_asym, quarantine, prob_mask, eff_mask, after_symptom, all);

		case_data = out.case_data;
		//total_cases = case_data.array_size;
		total_cases = out.cases_symptomatic;
		latest_onset = *max_element(case_data.onset.begin(), case_data.onset.end());
		extinct = *min_element(case_data.isolated.begin(), case_data.isolated.end());
		model.effective_R0.push_back(out.effective_R0);
		model.cases_in_gen.push_back(out.cases_in_gen);
		model.latest_onset.push_back(latest_onset);

		/*
		cout << "\nstep " << step << endl;
		//if (step == 1)
			show_data(out);
		/*else
		{
			cout << "effective R0 = " << out.effective_R0 << endl;
			cout << "number of cases in gen = " << out.cases_in_gen << endl;
		}
		cout << "total number of cases = " << total_cases << endl;
		cout << "extinct = " << extinct << endl;
		step++;
		cin.get();*/
	}

	model.mean_effective_R0 = 0.0;
	for (int i = 0; i < model.effective_R0.size(); i++)
		model.mean_effective_R0 = model.mean_effective_R0 + model.effective_R0[i];

	model.mean_effective_R0 = model.mean_effective_R0 / model.effective_R0.size();
	model.total_cases = total_cases;
	model.total_symptomatic = out.cases_symptomatic;
	model.total_asymptomatic = out.cases_asymptomatic;

	// find weekly cases
	int max_week = max_days / 7; // start to count from 0 --> 0 week, 1 week .. max_week week
	model.weekly_cases = vector<int>(max_week + 1, 0);
	int week = 0;
	for (int i = 0; i < case_data.array_size; i++)
	{
		if (!case_data.asym[i]) // if it is symptomatic
		{
			week = floor(case_data.onset[i] / 7);
			if (week <= max_week)
				model.weekly_cases[week] = model.weekly_cases[week] + 1;
		}

	}
	// cumulative number of cases
	int cumulative = 0;
	// define controlled outbreak
	int start_week = 12;
	int end_week = 16;
	vector<bool> control;
	model.controlled = false;
	for (int i = 0; i <= end_week; i++) // i = # week
	{
		cumulative = cumulative + model.weekly_cases[i];

		if (i >= start_week && i <= end_week)// 
		{
			if (model.weekly_cases[i] == 0 && cumulative < max_cases)
				control.push_back(true);
			else
				control.push_back(false);
		}
	}

	model.controlled = *min_element(control.begin(), control.end());

	return model;

}







int main()
{
	// Time since infectious
	double x_start = 0.0;	// start day
	double x_end = 15.0;	// end day
	double step = 0.01;		// step between days
	// Probability density function of distribution for parameters
	//delayDistribution(x_start, x_end, step);
	//incubationDistribution(x_start, x_end, step);
	//serialIntervalDistribution(x_start, x_end, step);

	/* Initial parameters */
	int num_sim = 1000;					// number of simulation to run
	int max_days = 365;					// max number of days for running process
	double num_initial_cases = 20;		// initial number of symptomatic cases
	double num_cases_asym = 20;			// initial number of asymptomatic cases (20)
	double max_cases = 100000.0;		// max number of cases for running process (5000.0 for outbreak control; 100000.0 for R_eff)
	double R0_isolated = 0.0;			// R0 for isolated cases
	double R0_community = 3.5;			// R0 for non-isolated cases
	double disp_iso = 1.0;				// dispersion parameter for negative binomial distribution for isolated cases
	double disp_com = 0.16;				// dispersion parameter for negative binomial distribution for non-isolated cases
	double k = 0.324;					// proportion of transmission before symptoms onset (1.95 -> 15%); (0.7 -> 30%); (30 -> <1%) (0.324 -> 40%)
	double delay_shape = 1.651524;		// shape of distribution for delay between symptom onset and isolation
	double delay_scale = 4.287786;		// scale of distribution for delay between symptom onset and isolation // 2.287786 --> mean=2 // 4.287786 --> mean = 3.8
	double prob_ascertain = 0.4;		// probability that cases are ascertained by contact tracing
	double prob_asym = 0.15;			// numerical proportion of cases that are subclinical (asymptotic): [0.0; 1.0] (0.15)
	double prob_asym_asym = 0.5;		// numerical proportion of asymptotic cases from asymptotic cases (0.5)
	double inf_asym = 0.6;				// infectivity of asymptotic case from initial R0 (0.6 -> 2.5 * 0.6 = 1.5)
	double Ra = R0_community * inf_asym;// reproduction number for asymptotic cases
	bool quarantine = false;			// quaratine is in effect or not (true --> traced contacts are isolated before symptoms onset)

	bool after_symptom = false;		// when people start to wear mask
	double prob_mask = 0.9;			// probability of wearing mask
	double eff_mask = 0.5;			// efficiency of mask wearing
	bool all = true;				// all people (infected and suspected) wear mask

	vector<Model_output> results;
	struct Model_output model;
	double mean_of_mean = 0.0;
	double mean_of_cases = 0.0;
	double outbreak_control = 0.0;


	ofstream fout_R0;
	//fout_R0.open("mean_eff_R0_tracrd80%.csv"); // 3)
	//fout_R0.open("masks_all70%_eff50%_R2.5_traced0%.csv");

	fout_R0.open("R3_sym20_asym20_prob70_tracing_40%.csv");
	// parameters
	fout_R0 << "num_sim" << "," << "max_days" << "," << "num_sym_cases" << "," << "num_asym_cases" << "," << "max_cases" << "," << "R0_isolated" << "," << "R0_community" << "," << "disp_iso" << ",";
	fout_R0 << "disp_com" << "," << "trans_before" << "," << "delay_shape" << "," << "delay_scale" << "," << "prob_ascertain" << "," << "prob_asym" << "," << "prob_asym_asym" << ",";
	fout_R0 << "Ra" << "," << "quarantine" << "," << "prob_mask" << "," << "eff_mask" << "," << "after_symptom" << "," << "all" << "\n";
	fout_R0 << num_sim << "," << max_days << "," << num_initial_cases << "," << num_cases_asym << "," << max_cases << "," << R0_isolated << "," << R0_community << "," << disp_iso << "," << disp_com << ",";
	fout_R0 << k << "," << delay_shape << "," << delay_scale << "," << prob_ascertain << "," << prob_asym << "," << prob_asym_asym << "," << Ra << "," << quarantine << ",";
	fout_R0 << prob_mask << "," << eff_mask << "," << after_symptom << "," << all << "\n";

	fout_R0 << "\n#Simulation" << "," << "Effective R0" << "," << "Total sympt" << "," << "Total asympt" << "," << "Latest onset" << "," << "Outbreak control" << "\n";


	for (int i = 1; i <= num_sim; i++)
	{
		printf("\nSimulation #%d\n", i);

		model = outbreak_model(num_initial_cases, num_cases_asym, max_days, max_cases, R0_isolated, R0_community,
			disp_iso, disp_com, k, delay_shape, delay_scale, prob_ascertain, prob_asym, prob_asym_asym, inf_asym, quarantine, prob_mask, eff_mask, after_symptom, all);
		results.push_back(model);
	}

	mean_of_mean = 0.0;
	mean_of_cases = 0.0;
	outbreak_control = 0.0;
	for (int i = 0; i < num_sim; i++)
	{
		mean_of_mean = mean_of_mean + results[i].mean_effective_R0;
		mean_of_cases = mean_of_cases + results[i].total_cases;
		outbreak_control = outbreak_control + results[i].controlled;
		fout_R0 << i + 1 << "," << results[i].mean_effective_R0 << "," << results[i].total_symptomatic << "," << results[i].total_asymptomatic << "," << results[i].latest_onset.back() << "," << results[i].controlled << "\n";

	}

	mean_of_mean = mean_of_mean / num_sim;
	mean_of_cases = mean_of_cases / num_sim;
	outbreak_control = 100.0 * outbreak_control / (num_sim * 1.0);
	//fout_R0 << "\nmean Ef. R0" << "," << mean_of_mean << "," << "mean cases" << "," << mean_of_cases << "," << "Outbreak control" << "," << outbreak_control << "\n";
	fout_R0 << "\nmean Ef. R0" << "," << setprecision(8) << mean_of_mean << ",";
	fout_R0 << "Outbreak control" << "," << setprecision(8) << outbreak_control << "\n";
	fout_R0.close();



	return 0;
}