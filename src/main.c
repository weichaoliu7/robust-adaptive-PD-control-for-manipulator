#include <stdio.h>
#include <stdlib.h>
#include "sine.h"
#include "cosine.h"

#define PI 3.14159

int main()
{

    double Ts = 0.0006; // sampling period
    double t0 = 0.0;
    double t1 = 3.0;
    double q1d_amp = 1.0, q2d_amp = 1.0;
    double dq1d_amp = 2.0 * PI, dq2d_amp = 2.0 * PI;
    double ddq1d_amp = -pow(2.0 * PI, 2), ddq2d_amp = -pow(2.0 * PI, 2);

    Data q1_d, q2_d, ddq1_d, ddq2_d;
    Data1 dq1_d, dq2_d;
    sine(&q1_d, Ts, t0, t1, q1d_amp); // target position of link 1
    sine(&q2_d, Ts, t0, t1, q2d_amp); // target position of link 2
    cosine(&dq1_d, Ts, t0, t1, dq1d_amp);
    cosine(&dq2_d, Ts, t0, t1, dq2d_amp);
    sine(&ddq1_d, Ts, t0, t1, ddq1d_amp);
    sine(&ddq2_d, Ts, t0, t1, ddq2d_amp);

    int ARRAY_SIZE = (t1 - t0) / Ts;

    // controller
    double ctrl_u1[ARRAY_SIZE], ctrl_u2[ARRAY_SIZE], ctrl_u3[ARRAY_SIZE]; // controller input
    double ctrl_u4[ARRAY_SIZE], ctrl_u5[ARRAY_SIZE], ctrl_u6[ARRAY_SIZE];
    double ctrl_u7[ARRAY_SIZE], ctrl_u8[ARRAY_SIZE], ctrl_u9[ARRAY_SIZE], ctrl_u10[ARRAY_SIZE];
    double x1[ARRAY_SIZE], x2[ARRAY_SIZE], x3[ARRAY_SIZE], x4[ARRAY_SIZE];

    ctrl_u1[0] = q1_d.y[0]; // system input
    ctrl_u2[0] = dq1_d.y[0];
    ctrl_u3[0] = ddq1_d.y[0];
    ctrl_u4[0] = q2_d.y[0];
    ctrl_u5[0] = dq2_d.y[0];
    ctrl_u6[0] = ddq2_d.y[0];

    ctrl_u7[0] = 0.1; // feedback from plant
    ctrl_u8[0] = 0;
    ctrl_u9[0] = 0.1;
    ctrl_u10[0] = 0;

    x1[0] = ctrl_u7[0]; // state variable
    x2[0] = ctrl_u8[0];
    x3[0] = ctrl_u9[0];
    x4[0] = ctrl_u10[0];

    double dq_d[2], ddq_d[2];
    double e[2], de[2];
    double y[2], dqr[2], ddqr[2];

    dq_d[0] = dq1_d.y[0];
    dq_d[1] = dq2_d.y[0];
    ddq_d[0] = ddq1_d.y[0];
    ddq_d[1] = ddq2_d.y[0];
    e[0] = ctrl_u7[0] - q1_d.y[0]; // position error
    e[1] = ctrl_u9[0] - q2_d.y[0];
    de[0] = ctrl_u8[0] - dq1_d.y[0]; // position error's derivative
    de[1] = ctrl_u10[0] - dq2_d.y[0];

    double gamma[2][2] = {{5, 0}, {0, 5}};
    for (int i = 0; i < 2; i++){
        y[i] = gamma[i][0] * e[0] + gamma[i][1] * e[1] + de[i]; // (2.11)
    }

    for (int i = 0; i < 2; i++){
        dqr[i] = dq_d[i] - gamma[i][0] * e[0] - gamma[i][1] * e[1]; // (2.12)
        ddqr[i] = ddq_d[i] - gamma[i][0] * de[0] - gamma[i][1] * de[1];
    }

    double d1 = 2.0, d2 = 3.0, d3 = 6.0;
    double e_norm, de_norm, sgn_y1, sgn_y2;
    e_norm = sqrt(e[0] * e[0] + e[1] * e[1]);
    de_norm = sqrt(de[0] * de[0] + de[1] * de[1]);
    sgn_y1 = y[0] >= 0 ? 1.0 : -1.0;
    sgn_y2 = y[1] >= 0 ? 1.0 : -1.0;

    double u1[ARRAY_SIZE], u2[ARRAY_SIZE];
    u1[0] = -(d1 + d2 * e_norm + d3 * de_norm) * sgn_y1; // (2.17)
    u2[0] = -(d1 + d2 * e_norm + d3 * de_norm) * sgn_y2;

    printf("u1[0] = %f\n", u1[0]);
    printf("u2[0] = %f\n", u2[0]);

    double Kp1[2][2] = {{80, 0}, {0, 90}};
    double Kp2[2][2] = {{50, 0}, {0, 50}};
    double Kv1[2][2] = {{80, 0}, {0, 80}};
    double Kv2[2][2] = {{50, 0}, {0, 50}};

    double alpha1 = 1.0, alpha2 = 1.0;
    double beta1 = 1.0, beta2 = 1.0;

    // adaptive regulator
    double g = 9.8, r1 = 1.0;
    double e1 = g / r1;

    double phi11, phi12, phi13, phi21, phi22, phi23;
    double phi[2][3], phi_T[3][3];
    double P[3][ARRAY_SIZE], dP[3], phi_T_y[3];
    double p1, p2, p3;

    p1 = 4.1;
    p2 = 1.9;
    p3 = 1.7;

    P[0][0] = p1;
    P[1][0] = p2;
    P[2][0] = p3;

    phi11 = ddqr[0] + e1 * cos(x3[0]); // regression matrix
    phi12 = ddqr[0] + ddqr[1];
    phi13 = 2 * ddqr[0] * cos(x3[0]) + ddqr[1] * cos(x3[0]) - x4[0] * dqr[0] * sin(x3[0]) - (x2[0] + x4[0]) * dqr[1] * sin(x3[0]) + e1 * cos(x1[0] + x3[0]);
    phi21 = 0;
    phi22 = phi12;
    phi23 = x2[0] * dqr[0] * sin(x3[0]) + ddqr[0] * cos(x3[0]) + e1 * cos(x1[0] + x3[0]);

    phi_T[0][0] = phi11;
    phi_T[1][0] = phi12;
    phi_T[2][0] = phi13;
    phi_T[0][1] = phi21;
    phi_T[1][1] = phi22;
    phi_T[2][1] = phi23;
    phi_T[0][2] = phi_T[1][2] = phi_T[3][2] = 0;

    phi[0][0] = phi11;
    phi[0][1] = phi12;
    phi[0][2] = phi13;
    phi[1][0] = phi21;
    phi[1][1] = phi22;
    phi[1][2] = phi23;

    double Gamma[3][3] = {{5, 0, 0}, {0, 5, 0}, {0, 0, 5}};

    double Gamma_phi_T[3][3];
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            double sum = 0;
            for (int k = 0; k < 3; k++){
                sum += Gamma[i][k] * phi_T[k][j];
            }
            Gamma_phi_T[i][j] = sum;
        }
    }

    for (int i = 0; i < 3; i++){
        dP[i] = - (Gamma_phi_T[i][0] * y[0] + Gamma_phi_T[i][1] * y[1]);
    }

    p1 = P[0][0] + dP[0] * Ts;
    p2 = P[1][0] + dP[1] * Ts;
    p3 = P[2][0] + dP[2] * Ts;

    P[0][0] = p1;
    P[1][0] = p2;
    P[2][0] = p3;

    double phi_P[2];

    for (int i = 0; i < 2; i++) {
        phi_P[i] = 0;
        for (int j = 0; j < 3; j++) {
            phi_P[i] += phi[i][j] * P[j][0];
        }
    }

    double tol[2][ARRAY_SIZE];
    double tol1[ARRAY_SIZE], tol2[ARRAY_SIZE];

    tol[0][0] = -(Kp1[0][0] + Kp2[0][0] / (alpha1 + fabs(e[0]))) * e[0] -
                (Kv1[0][0] + Kv2[0][0] / (beta1 + fabs(de[0]))) * de[0] +
                phi_P[0] + u1[0]; // (2.16)

    tol[1][0] = -(Kp1[1][1] + Kp2[1][1] / (alpha2 + fabs(e[1]))) * e[1] -
                (Kv1[1][1] + Kv2[1][1] / (beta2 + fabs(de[1]))) * de[1] +
                phi_P[1] + u2[0];

    tol1[0] = tol[0][0];
    tol2[0] = tol[1][0];

    double ctrl_y1[ARRAY_SIZE], ctrl_y2[ARRAY_SIZE];
    ctrl_y1[0] = tol1[0];
    ctrl_y2[0] = tol2[0];

    // plant
    double plant_u1[ARRAY_SIZE], plant_u2[ARRAY_SIZE];
    plant_u1[0] = ctrl_y1[0];
    plant_u2[0] = ctrl_y2[0];

    double m1 = 0.5, m2 = 0.5, r2 = 0.8;
    double D11, D12, D22, C12;
    double D[2][2], C[2][2];

    D11 = (m1 + m2) * pow(r1, 2) + m2 * pow(r2, 2) + 2 * m2 * r1 * r2 * cos(x3[0]);
    D12 = m2 * pow(r2, 2) + m2 * r1 * r2 * cos(x3[0]);
    D22 = m2 * pow(r2, 2);
    D[0][0] = D11;
    D[0][1] = D12;
    D[1][0] = D12;
    D[1][1] = D22;

    C12 = m2 * r1 * r2 * sin(x3[0]);
    C[0][0] = -C12 * x4[0];
    C[0][1] = -C12 * (x2[0] + x4[0]);
    C[1][0] = C12 * x2[0];
    C[1][1] = 0;

    double g1, g2, G[2], w1, w2, w[2];
    g1 = (m1 + m2) * r1 * cos(x3[0]) + m2 * r2 * cos(x1[0] + x3[0]);
    g2 = m2 * r2 * cos(x1[0] + x3[0]);

    G[0] = g1 * g;
    G[1] = g2 * g;

    w1 = 1.5 + 2.0 * e[0] + 5.0 * de[0];
    w2 = 1.5 + 2.0 * e[1] + 5.0 * de[1];
    w[0] = w1;
    w[1] = w2;

    double a = D[0][0], b = D[0][1], c = D[1][0], d = D[1][1];
    double D_inv[2][2], S[2];

    D_inv[0][0] = d / (a * d - b * c);
    D_inv[0][1] = -b / (a * d - b * c);
    D_inv[1][0] = -c / (a * d - b * c);
    D_inv[1][1] = a / (a * d - b * c);

    double dq[2][ARRAY_SIZE], ddq[2][ARRAY_SIZE], q[2][ARRAY_SIZE];
    q[0][0] = x1[0];
    q[1][0] = x3[0];
    dq[0][0] = x2[0];
    dq[1][0] = x4[0];

    double tol_Cdq_G_w[2], tol_Cdq_G_w1, tol_Cdq_G_w2;
    tol_Cdq_G_w1 = tol[0][0] - (C[0][0] * dq[0][0] + C[0][1] * dq[1][0]) - G[0] - w[0];
    tol_Cdq_G_w2 = tol[1][0] - (C[1][0] * dq[0][0] + C[1][1] * dq[1][0]) - G[1] - w[1];

    S[0] = D_inv[0][0] * tol_Cdq_G_w1 + D_inv[0][1] * tol_Cdq_G_w2;
    S[1] = D_inv[1][0] * tol_Cdq_G_w1 + D_inv[1][1] * tol_Cdq_G_w2;

    ddq[0][0] = S[0]; // second order position derivative of link 1
    ddq[1][0] = S[1]; // second order position derivative of link 2

    dq[0][0] = dq[0][0] + ddq[0][0] * Ts;
    dq[1][0] = dq[1][0] + ddq[1][0] * Ts;

    q[0][0] = q[0][0] + dq[0][0] * Ts;
    q[1][0] = q[1][0] + dq[1][0] * Ts;

    double plant_y1[ARRAY_SIZE], plant_y2[ARRAY_SIZE], plant_y3[ARRAY_SIZE], plant_y4[ARRAY_SIZE];

    plant_y1[0] = q[0][0];
    plant_y2[0] = dq[0][0];
    plant_y3[0] = q[1][0];
    plant_y4[0] = dq[1][0];


    for (int i = 1; i < ARRAY_SIZE; i++)
    {

        // controller
        ctrl_u1[i] = q1_d.y[i]; // system input
        ctrl_u2[i] = dq1_d.y[i];
        ctrl_u3[i] = ddq1_d.y[i];
        ctrl_u4[i] = q2_d.y[i];
        ctrl_u5[i] = dq2_d.y[i];
        ctrl_u6[i] = ddq2_d.y[i];

        ctrl_u7[i] = plant_y1[i -1]; // feedback from plant
        ctrl_u8[i] = plant_y2[i -1];
        ctrl_u9[i] = plant_y3[i -1];
        ctrl_u10[i] = plant_y4[i -1];

        x1[i] = ctrl_u7[i]; // state variable
        x2[i] = ctrl_u8[i];
        x3[i] = ctrl_u9[i];
        x4[i] = ctrl_u10[i];

        dq_d[0] = dq1_d.y[i];
        dq_d[1] = dq2_d.y[i];
        ddq_d[0] = ddq1_d.y[i];
        ddq_d[1] = ddq2_d.y[i];
        e[0] = ctrl_u7[i] - q1_d.y[i]; // position error
        e[1] = ctrl_u9[i] - q2_d.y[i];
        de[0] = ctrl_u8[i] - dq1_d.y[i]; // position error's derivative
        de[1] = ctrl_u10[i] - dq2_d.y[i];

        for (int j = 0; j < 2; j++){
            y[j] = gamma[j][0] * e[0] + gamma[j][1] * e[1] + de[j]; // (2.11)
        }

        for (int j = 0; j < 2; j++){
            dqr[j] = dq_d[j] - gamma[j][0] * e[0] - gamma[j][1] * e[1]; // (2.12)
            ddqr[j] = ddq_d[j] - gamma[j][0] * de[0] - gamma[j][1] * de[1];
        }

        e_norm = sqrt(e[0] * e[0] + e[1] * e[1]);
        de_norm = sqrt(de[0] * de[0] + de[1] * de[1]);
        sgn_y1 = y[0] >= 0 ? 1.0 : -1.0;
        sgn_y2 = y[1] >= 0 ? 1.0 : -1.0;

        u1[i] = -(d1 + d2 * e_norm + d3 * de_norm) * sgn_y1; // (2.17)
        u2[i] = -(d1 + d2 * e_norm + d3 * de_norm) * sgn_y2;

        // adaptive regulator
        phi11 = ddqr[0] + e1 * cos(x3[i]); // regression matrix
        phi12 = ddqr[0] + ddqr[1];
        phi13 = 2 * ddqr[0] * cos(x3[i]) + ddqr[1] * cos(x3[i]) - x4[i] * dqr[0] * sin(x3[i]) - (x2[i] + x4[i]) * dqr[1] * sin(x3[i]) + e1 * cos(x1[i] + x3[i]);
        phi21 = 0;
        phi22 = phi12;
        phi23 = x2[i] * dqr[0] * sin(x3[i]) + ddqr[0] * cos(x3[i]) + e1 * cos(x1[i] + x3[i]);

        phi_T[0][0] = phi11;
        phi_T[1][0] = phi12;
        phi_T[2][0] = phi13;
        phi_T[0][1] = phi21;
        phi_T[1][1] = phi22;
        phi_T[2][1] = phi23;
        phi_T[0][2] = phi_T[1][2] = phi_T[3][2] = 0;

        phi[0][0] = phi11;
        phi[0][1] = phi12;
        phi[0][2] = phi13;
        phi[1][0] = phi21;
        phi[1][1] = phi22;
        phi[1][2] = phi23;

        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                double sum = 0;
                for (int l = 0; l < 3; l++){
                    sum += Gamma[j][l] * phi_T[l][k];
                }
                Gamma_phi_T[j][k] = sum;
            }
        }

        for (int j = 0; j < 3; j++){
            dP[j] = - (Gamma_phi_T[j][0] * y[0] + Gamma_phi_T[j][1] * y[1]);
        }

        p1 = P[0][i - 1] + dP[0] * Ts;
        p2 = P[1][i - 1] + dP[1] * Ts;
        p3 = P[2][i - 1] + dP[2] * Ts;

        P[0][i] = p1;
        P[1][i] = p2;
        P[2][i] = p3;

        for (int j = 0; j < 2; j++) {
            phi_P[j] = 0;
            for (int k = 0; k < 3; k++) {
                phi_P[j] += phi[j][k] * P[k][i];
            }
        }

        tol[0][i] = -(Kp1[0][0] + Kp2[0][0] / (alpha1 + fabs(e[0]))) * e[0] -
                    (Kv1[0][0] + Kv2[0][0] / (beta1 + fabs(de[0]))) * de[0] +
                    phi_P[0] + u1[i]; // (2.16)

        tol[1][i] = -(Kp1[1][1] + Kp2[1][1] / (alpha2 + fabs(e[1]))) * e[1] -
                    (Kv1[1][1] + Kv2[1][1] / (beta2 + fabs(de[1]))) * de[1] +
                    phi_P[1] + u2[i];

        tol1[i] = tol[0][i];
        tol2[i] = tol[1][i];

        ctrl_y1[i] = tol1[i];
        ctrl_y2[i] = tol2[i];

        // plant
        plant_u1[i] = ctrl_y1[i];
        plant_u2[i] = ctrl_y2[i];

        D11 = (m1 + m2) * pow(r1, 2) + m2 * pow(r2, 2) + 2 * m2 * r1 * r2 * cos(x3[i]);
        D12 = m2 * pow(r2, 2) + m2 * r1 * r2 * cos(x3[i]);
        D22 = m2 * pow(r2, 2);
        D[0][0] = D11;
        D[0][1] = D12;
        D[1][0] = D12;
        D[1][1] = D22;

        C12 = m2 * r1 * r2 * sin(x3[i]);
        C[0][0] = -C12 * x4[i];
        C[0][1] = -C12 * (x2[i] + x4[i]);
        C[1][0] = C12 * x2[0];
        C[1][1] = 0;

        g1 = (m1 + m2) * r1 * cos(x3[i]) + m2 * r2 * cos(x1[i] + x3[i]);
        g2 = m2 * r2 * cos(x1[i] + x3[i]);

        G[0] = g1 * g;
        G[1] = g2 * g;

        w1 = 1.5 + 2.0 * e[0] + 5.0 * de[0];
        w2 = 1.5 + 2.0 * e[1] + 5.0 * de[1];
        w[0] = w1;
        w[1] = w2;

        a = D[0][0], b = D[0][1], c = D[1][0], d = D[1][1];

        D_inv[0][0] = d / (a * d - b * c);
        D_inv[0][1] = -b / (a * d - b * c);
        D_inv[1][0] = -c / (a * d - b * c);
        D_inv[1][1] = a / (a * d - b * c);

        q[0][i] = x1[i];
        q[1][i] = x3[i];
        dq[0][i] = x2[i];
        dq[1][i] = x4[i];

        tol_Cdq_G_w1 = tol[0][i] - (C[0][0] * dq[0][i] + C[0][1] * dq[1][i]) - G[0] - w[0];
        tol_Cdq_G_w2 = tol[1][i] - (C[1][0] * dq[0][i] + C[1][1] * dq[1][i]) - G[1] - w[1];

        S[0] = D_inv[0][0] * tol_Cdq_G_w1 + D_inv[0][1] * tol_Cdq_G_w2;
        S[1] = D_inv[1][0] * tol_Cdq_G_w1 + D_inv[1][1] * tol_Cdq_G_w2;

        ddq[0][i] = S[0]; // second order position derivative of link 1
        ddq[1][i] = S[1]; // second order position derivative of link 2

        dq[0][i] = dq[0][i - 1] + ddq[0][i] * Ts;
        dq[1][i] = dq[1][i - 1] + ddq[1][i] * Ts;

        q[0][i] = q[0][i - 1] + dq[0][i] * Ts;
        q[1][i] = q[1][i - 1] + dq[1][i] * Ts;

        plant_y1[i] = q[0][i];
        plant_y2[i] = dq[0][i];
        plant_y3[i] = q[1][i];
        plant_y4[i] = dq[1][i];

    }

    FILE *fp1 = fopen("q1.txt", "w");
    if (fp1 == NULL){
        printf("Failed to open file q1.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp1, "%.6f\n", plant_y1[i]);
    }
    fclose(fp1);
    printf("q1.txt successfully saved\n");

    FILE *fp2 = fopen("q2.txt", "w");
    if (fp2 == NULL){
        printf("Failed to open file q2.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp1, "%.6f\n", plant_y3[i]);
    }
    fclose(fp2);
    printf("q2.txt successfully saved\n");

    FILE *fp3 = fopen("dq1.txt", "w");
    if (fp3 == NULL){
        printf("Failed to open file dq1.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp3, "%.6f\n", plant_y2[i]);
    }
    fclose(fp3);
    printf("dq1.txt successfully saved\n");

    FILE *fp4 = fopen("dq2.txt", "w");
    if (fp4 == NULL){
        printf("Failed to open file dq2.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp4, "%.6f\n", plant_y4[i]);
    }
    fclose(fp4);
    printf("dq2.txt successfully saved\n");

    FILE *fp5 = fopen("tol1.txt", "w");
    if (fp5 == NULL){
        printf("Failed to open file tol1.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp5, "%.6f\n", tol[0][i]);
    }
    fclose(fp5);
    printf("tol1.txt successfully saved\n");

    FILE *fp6 = fopen("tol2.txt", "w");
    if (fp6 == NULL){
        printf("Failed to open file tol2.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp6, "%.6f\n", tol[1][i]);
    }
    fclose(fp6);
    printf("to12.txt successfully saved\n");

    FILE *fp7 = fopen("p1.txt", "w");
    if (fp7 == NULL){
        printf("Failed to open file p1.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp7, "%.6f\n", P[0][i]);
    }
    fclose(fp7);
    printf("p1.txt successfully saved\n");

    FILE *fp8 = fopen("p2.txt", "w");
    if (fp8 == NULL){
        printf("Failed to open file p2.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp8, "%.6f\n", P[1][i]);
    }
    fclose(fp8);
    printf("p2.txt successfully saved\n");

    FILE *fp9 = fopen("p3.txt", "w");
    if (fp9 == NULL){
        printf("Failed to open file p3.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++){
        fprintf(fp9, "%.6f\n", P[2][i]);
    }
    fclose(fp9);
    printf("p3.txt successfully saved\n");

    // set title "plot of q1"
    // plot "q1.txt" with lines

    return 0;
}
