/*
	Generalized reconstruction algorithm developed for implementation on ZedBoard.
	Inputs from MATLAB: measurement matrix phi, dictionaries for digits '0' through '9'
	Input from the board: measured signal
	Output: reconstructed data, reconstruction coefficients
	Written by: Karlo Sever
*/

#include "measurement_matrix.h"
#include "eigv_matrix_0.h"
#include "eigv_matrix_1.h"
#include "eigv_matrix_2.h"
#include "eigv_matrix_3.h"
#include "eigv_matrix_4.h"
#include "eigv_matrix_5.h"
#include "eigv_matrix_6.h"
#include "eigv_matrix_7.h"
#include "eigv_matrix_8.h"
#include "eigv_matrix_9.h"
#include "test6.h"

#define DIM 784													//Image dimensions 28 x 28 = 784
#define MDIM 78													//Length of measured signal
#define NDIM 10													//Number of digits

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main(void) {

	VectorXd x_test(DIM);
	for (int i = 0; i < DIM; i++)
		x_test(i) = xTest6[i];									//Input test image to Eigen library VectorXd

	MatrixXd phi(MDIM, DIM);
	int indx = 0;
	std::cout << "Acquiring measurement matrix phi... \n\n\n";
	for (int i = 0; i < MDIM; i++) {
		for (int j = 0; j < DIM; j++) {
			phi(i, j) = phi_vect[indx++];						//Input measurement Eigen library MatrixXd
		}
	}

	MatrixXd eigv0(DIM, DIM), eigv1(DIM, DIM), eigv2(DIM, DIM), eigv3(DIM, DIM), eigv4(DIM, DIM), 
			 eigv5(DIM, DIM), eigv6(DIM, DIM), eigv7(DIM, DIM), eigv8(DIM, DIM), eigv9(DIM, DIM);
	indx = 0;
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {							//Input dictionary to Eigen library MatrixXd
			eigv0(i, j) = v0[indx];
			eigv1(i, j) = v1[indx];
			eigv2(i, j) = v2[indx];
			eigv3(i, j) = v3[indx];
			eigv4(i, j) = v4[indx];
			eigv5(i, j) = v5[indx];
			eigv6(i, j) = v6[indx];
			eigv7(i, j) = v7[indx];
			eigv8(i, j) = v8[indx];
			eigv9(i, j) = v9[indx];
			indx++;
		}
	}

	for (int i = 0; i < DIM; i++) {
		printf("%6.2f", x_test[i]);								//Check original image (visual output)
		if (i % 28 == 0)
			printf("\n");
	}
	
	VectorXd y(MDIM);
	y = phi * x_test;
	std::cout << "\nMeasured signal y... \n" << y << std::endl;

																//Reconstruction based on compressed sensing
	VectorXd alpha_rec(DIM), x_rec(DIM), recon(DIM);
	double arec[DIM][NDIM], xrec[DIM][NDIM];
	double arecL1[10], minL1;
	int minL1_digit = 0;

	x_rec = phi.transpose()*(phi*x_test);

	for (int i = 0; i < NDIM; i++) {

		switch (i) {
			case 0: alpha_rec = eigv0.transpose() * x_rec;
					break;
			case 1: alpha_rec = eigv1.transpose() * x_rec;
					break;
			case 2: alpha_rec = eigv2.transpose() * x_rec;
					break;
			case 3: alpha_rec = eigv3.transpose() * x_rec;
					break;
			case 4: alpha_rec = eigv4.transpose() * x_rec;
					break;
			case 5: alpha_rec = eigv5.transpose() * x_rec;
					break;
			case 6: alpha_rec = eigv6.transpose() * x_rec;
					break;
			case 7: alpha_rec = eigv7.transpose() * x_rec;
					break;
			case 8: alpha_rec = eigv8.transpose() * x_rec;
					break;
			case 9: alpha_rec = eigv9.transpose() * x_rec;
					break;
		}
		
		for (int k = 0; k < DIM; k++) {
			arec[k][i] = alpha_rec(k);
		}

		arecL1[i] = alpha_rec.template lpNorm<1>();				//Eigen library provides a way of computing L1-norm
		std::cout << "\n\nL1-norm of coefficients of representation of the signal over different dictionaries\n" << arecL1[i];
		
	}
	
	minL1 = arecL1[0];
	for (int i = 1; i < NDIM; i++) {
		if (arecL1[i] < minL1) {
			minL1 = arecL1[i];
			minL1_digit = i;
		}
	}
	std::cout << "\n\nMinimal L1-norm is achieved with the representation of the original image with the dictionary of digit " << minL1_digit << ".\n\n";

	for (int k = 0; k < DIM; k++) {
		alpha_rec(k) = arec[k][minL1_digit];
	}

	switch (minL1_digit) {
		case 0: recon = eigv0 * alpha_rec;
			break;
		case 1: recon = eigv1 * alpha_rec;
			break;
		case 2: recon = eigv2 * alpha_rec;
			break;
		case 3: recon = eigv3 * alpha_rec;
			break;
		case 4: recon = eigv4 * alpha_rec;
			break;
		case 5: recon = eigv5 * alpha_rec;
			break;
		case 6: recon = eigv6 * alpha_rec;
			break;
		case 7: recon = eigv7 * alpha_rec;
			break;
		case 8: recon = eigv8 * alpha_rec;
			break;
		case 9: recon = eigv9 * alpha_rec;
			break;
	}

	for (int i = 0; i < DIM; i++) {
		printf("%6.2f ", recon[i]);
		if (i % 28 == 0) {
			printf("\n");
		}
	}
	std::cout << std::endl;

	return 0;
}
