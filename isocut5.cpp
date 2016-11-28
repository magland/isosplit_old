#include "isocut5.h"
#include "jisotonic5.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

namespace ns_isocut5 {
void copy_samples(long N, float* out, float* in);
double sum(long N, float* X);
long find_min_index(long N, float* X);
long find_max_index(long N, float* X);
double compute_ks4(long N, float* counts1, float* counts2);
double compute_ks5(long* critical_range_min, long* critical_range_max, long N, float* counts1, float* counts2, long peak_index);
void debug_print_array(long N, float* X);
}

void isocut5_mex(double* dipscore, double* cutpoint, int N, double* samples)
{
    *dipscore = 0;
    *cutpoint = 0;
    isocut5_opts opts;
    float* samplesf = (float*)malloc(sizeof(float) * N);
    for (long i = 0; i < N; i++)
        samplesf[i] = samples[i];
    isocut5(dipscore, cutpoint, N, samplesf, opts);
    free(samplesf);
}

void isocut5(double* dipscore_out, double* cutpoint_out, long N, float* samples, isocut5_opts opts)
{
    float* samples_sorted = (float*)malloc(sizeof(float) * N);

    // sort the samples if needed
    if (opts.already_sorted)
        ns_isocut5::copy_samples(N, samples_sorted, samples);
    else
        jisotonic5_sort(N, samples_sorted, samples);

    float num_bins_factor = 1;
    long num_bins = ceil(sqrt(N * 1.0 / 2) * num_bins_factor);

    long num_bins_1 = ceil(num_bins / 2);
    long num_bins_2 = num_bins - num_bins_1;
    long num_intervals = num_bins_1 + num_bins_2;
    float* intervals = (float*)malloc(sizeof(float) * num_intervals);
    for (int i = 0; i < num_bins_1; i++)
        intervals[i] = i + 1;
    for (int i = 0; i < num_bins_2; i++)
        intervals[num_intervals - 1 - i] = i + 1;
    float alpha = (N - 1) / ns_isocut5::sum(num_intervals, intervals);
    for (int i = 0; i < num_intervals; i++)
        intervals[i] *= alpha;
    long N_sub = num_intervals + 1;
    float* inds = (float*)malloc(sizeof(float) * N_sub);
    inds[0] = 0;
    for (long i = 0; i < num_intervals; i++)
        inds[i + 1] = inds[i] + intervals[i];
    float* X_sub = (float*)malloc(sizeof(float) * N_sub);
    for (long i = 0; i < N_sub; i++)
        X_sub[i] = samples_sorted[(long)inds[i]];
    float* densities = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* spacings = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* multiplicities = (float*)malloc(sizeof(float) * (N_sub - 1));
    for (long i = 0; i < N_sub - 1; i++) {
        spacings[i] = X_sub[i + 1] - X_sub[i];
        multiplicities[i] = ((long)inds[i + 1]) - ((long)inds[i]);
        densities[i] = multiplicities[i] / spacings[i];
    }

    float* densities_unimodal_fit = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* densities_resid = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* densities_unimodal_fit_times_spacings = (float*)malloc(sizeof(float) * (N_sub - 1));
    jisotonic5_updown(N_sub - 1, densities_unimodal_fit, densities, multiplicities);
    for (long i = 0; i < N_sub - 1; i++)
        densities_resid[i] = densities[i] - densities_unimodal_fit[i];
    for (long i = 0; i < N_sub - 1; i++)
        densities_unimodal_fit_times_spacings[i] = densities_unimodal_fit[i] * spacings[i];
    long critical_range_min, critical_range_max;
    long peak_index = ns_isocut5::find_max_index(N_sub - 1, densities_unimodal_fit);
    *dipscore_out = ns_isocut5::compute_ks5(&critical_range_min, &critical_range_max, N_sub - 1, multiplicities, densities_unimodal_fit_times_spacings, peak_index);
    long critical_range_length = critical_range_max - critical_range_min + 1;

    float* densities_resid_on_critical_range = (float*)malloc(sizeof(float) * (critical_range_length));
    float* densities_resid_fit_on_critical_range = (float*)malloc(sizeof(float) * (critical_range_length));
    float* weights_for_downup = (float*)malloc(sizeof(float) * (critical_range_length));
    for (long i = 0; i < critical_range_length; i++) {
        densities_resid_on_critical_range[i] = densities_resid[critical_range_min + i];
        weights_for_downup[i] = spacings[critical_range_min + i];
    }
    jisotonic5_downup(critical_range_length, densities_resid_fit_on_critical_range, densities_resid_on_critical_range, weights_for_downup);

    long cutpoint_index = ns_isocut5::find_min_index(critical_range_length, densities_resid_fit_on_critical_range);
    *cutpoint_out = (X_sub[critical_range_min + cutpoint_index] + X_sub[critical_range_min + cutpoint_index + 1]) / 2;

    free(samples_sorted);
    free(densities_unimodal_fit);
    free(densities_resid);
    free(weights_for_downup);
    free(intervals);
    free(inds);
    free(X_sub);
    free(multiplicities);
    free(spacings);
    free(densities);
    free(densities_unimodal_fit_times_spacings);
    free(densities_resid_on_critical_range);
    free(densities_resid_fit_on_critical_range);
}

void isocut5_old(double* dipscore_out, double* cutpoint_out, long N, float* samples, isocut5_opts opts)
{
    float* samples_sorted = (float*)malloc(sizeof(float) * N);

    // sort the samples if needed
    if (opts.already_sorted)
        ns_isocut5::copy_samples(N, samples_sorted, samples);
    else
        jisotonic5_sort(N, samples_sorted, samples);

    float num_bins_factor = 1;
    long num_bins = ceil(sqrt(N * 1.0 / 2) * num_bins_factor);

    long num_bins_1 = ceil(num_bins / 2);
    long num_bins_2 = num_bins - num_bins_1;
    long num_intervals = num_bins_1 + num_bins_2;
    float* intervals = (float*)malloc(sizeof(float) * num_intervals);
    for (int i = 0; i < num_bins_1; i++)
        intervals[i] = i + 1;
    for (int i = 0; i < num_bins_2; i++)
        intervals[num_intervals - 1 - i] = i + 1;
    float alpha = (N - 1) / ns_isocut5::sum(num_intervals, intervals);
    for (int i = 0; i < num_intervals; i++)
        intervals[i] *= alpha;
    long N_sub = num_intervals + 1;
    float* inds = (float*)malloc(sizeof(float) * N_sub);
    inds[0] = 0;
    for (long i = 0; i < num_intervals; i++)
        inds[i + 1] = inds[i] + intervals[i];
    float* X_sub = (float*)malloc(sizeof(float) * N_sub);
    for (long i = 0; i < N_sub; i++)
        X_sub[i] = samples_sorted[(long)inds[i]];
    float* densities = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* spacings = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* multiplicities = (float*)malloc(sizeof(float) * (N_sub - 1));
    for (long i = 0; i < N_sub - 1; i++) {
        spacings[i] = X_sub[i + 1] - X_sub[i];
        multiplicities[i] = ((long)inds[i + 1]) - ((long)inds[i]);
        densities[i] = multiplicities[i] / spacings[i];
    }

    float* densities_unimodal_fit = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* densities_resid = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* densities_resid_fit = (float*)malloc(sizeof(float) * (N_sub - 1));
    float* weights_for_downup = (float*)malloc(sizeof(float) * (N_sub - 1));

    jisotonic5_updown(N_sub - 1, densities_unimodal_fit, densities, multiplicities);
    for (long i = 0; i < N_sub - 1; i++)
        densities_resid[i] = densities[i] - densities_unimodal_fit[i];
    for (long i = 0; i < N_sub - 1; i++)
        weights_for_downup[i] = spacings[i];
    jisotonic5_downup(N_sub - 1, densities_resid_fit, densities_resid, weights_for_downup);

    long cutpoint_index = ns_isocut5::find_min_index(N_sub - 1, densities_resid_fit);
    *cutpoint_out = (X_sub[cutpoint_index] + X_sub[cutpoint_index + 1]) / 2;

    float* densities_unimodal_fit_times_spacings = (float*)malloc(sizeof(float) * (N_sub - 1));
    for (long i = 0; i < N_sub - 1; i++)
        densities_unimodal_fit_times_spacings[i] = densities_unimodal_fit[i] * spacings[i];
    *dipscore_out = ns_isocut5::compute_ks4(N_sub - 1, multiplicities, densities_unimodal_fit_times_spacings);

    free(samples_sorted);
    free(densities_unimodal_fit);
    free(densities_resid);
    free(densities_resid_fit);
    free(weights_for_downup);
    free(intervals);
    free(inds);
    free(X_sub);
    free(multiplicities);
    free(spacings);
    free(densities);
    free(densities_unimodal_fit_times_spacings);
}

namespace ns_isocut5 {

void copy_samples(long N, float* out, float* in)
{
    for (long i = 0; i < N; i++)
        out[i] = in[i];
}

double sum(long N, float* X)
{
    double ret = 0;
    for (long i = 0; i < N; i++)
        ret += X[i];
    return ret;
}

long find_min_index(long N, float* X)
{
    long ret = 0;
    for (long i = 0; i < N; i++) {
        if (X[i] < X[ret])
            ret = i;
    }
    return ret;
}

long find_max_index(long N, float* X)
{
    long ret = 0;
    for (long i = 0; i < N; i++) {
        if (X[i] > X[ret])
            ret = i;
    }
    return ret;
}

double compute_ks4(long N, float* counts1, float* counts2)
{
    double sum_counts1 = sum(N, counts1);
    double sum_counts2 = sum(N, counts2);

    double cumsum_counts1 = 0;
    double cumsum_counts2 = 0;

    double max_diff = 0;
    for (long i = 0; i < N; i++) {
        cumsum_counts1 += counts1[i];
        cumsum_counts2 += counts2[i];
        double diff = fabs(cumsum_counts1 / sum_counts1 - cumsum_counts2 / sum_counts2);
        if (diff > max_diff)
            max_diff = diff;
    }

    return max_diff * sqrt((sum_counts1 + sum_counts2) / 2);
}

double compute_ks5(long* critical_range_min, long* critical_range_max, long N, float* counts1, float* counts2, long peak_index)
{
    *critical_range_min = 0;
    *critical_range_max = N - 1; //should get over-written!
    double ks_best = -1;

    // from the left
    {
        float* counts1_left = (float*)malloc(sizeof(float) * (peak_index + 1));
        float* counts2_left = (float*)malloc(sizeof(float) * (peak_index + 1));
        for (long i = 0; i <= peak_index; i++) {
            counts1_left[i] = counts1[i];
            counts2_left[i] = counts2[i];
        }
        long len = peak_index + 1;
        while ((len >= 4) || (len == peak_index + 1)) {
            double ks0 = compute_ks4(len, counts1_left, counts2_left);
            if (ks0 > ks_best) {
                *critical_range_min = 0;
                *critical_range_max = len - 1;
                ks_best = ks0;
            }
            len = len / 2;
        }
    }

    // from the right
    {
        float* counts1_right = (float*)malloc(sizeof(float) * (N - peak_index));
        float* counts2_right = (float*)malloc(sizeof(float) * (N - peak_index));
        for (long i = 0; i < N - peak_index; i++) {
            counts1_right[i] = counts1[N - 1 - i];
            counts2_right[i] = counts2[N - 1 - i];
        }
        long len = N - peak_index;
        while ((len >= 4) || (len == N - peak_index)) {
            double ks0 = compute_ks4(len, counts1_right, counts2_right);
            if (ks0 > ks_best) {
                *critical_range_min = N - len;
                *critical_range_max = N - 1;
                ks_best = ks0;
            }
            len = len / 2;
        }
    }

    return ks_best;
}

void debug_print_array(long N, float* X)
{
    for (long i = 0; i < N; i++) {
        if ((i > 0) && (i % 10 == 0))
            printf("\n");
        printf("%g ", X[i]);
    }
    printf("\n");
}
}
