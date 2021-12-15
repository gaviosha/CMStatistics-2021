# Robust inference for change points in piecewise polynomials of general degrees

This repository contains code for building the presentation and reproducing numerical examples in my talk at CMStatistics 2021. More information about the talk can be found [here](http://www.cmstatistics.org/RegistrationsV2/CMStatistics2021/viewSubmission.php?in=266&token=pp7561qnq930q8012s476r34o7r6rs66).

## Abstract

Multiple change-point detection has become popular with the routine collection of complex non-stationary time series. An equally important but comparatively neglected question concerns quantifying the level of uncertainty around each putative changepoint. Though a handful of procedures exist in the literature, most all make assumptions on the density of the contaminating noise which are impossible to verify in practice. We present a procedure that, under minimal assumptions, returns localized regions of a data sequence that must contain a changepoint at some global significance level chosen by the user. Our procedure is computationally efficient, applicable to change points in higher-order polynomials, and moreover, all results are fully non-asymptotic. We will discuss some appealing theoretical properties of our procedure, and show its good practical performance on real and simulated data.
