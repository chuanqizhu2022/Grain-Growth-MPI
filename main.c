// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/2/16
// Updated on 2022/06/20

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

#define NDX 100 //差分計算における計算領域一辺の分割数
#define NDY 100 //差分計算における計算領域一辺の分割数
#define NDZ 100
#define N 10 //考慮する結晶方位の数＋１(MPF0.cppと比較して、この値を大きくしている)
#define BEGIN 1
#define UTAG 2
#define DTAG 3
#define NONE 0
#define DONE 4
#define MASTER 0

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1; //計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
int ndmz = NDZ - 1;
int nm = N - 1, nmm = N - 2; //考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
double PI = 3.141592;        //π、計算カウント数
double RR = 8.3145;          //ガス定数

double aij[N][N]; //勾配エネルギー係数
double wij[N][N]; //ペナルティー項の係数
double mij[N][N]; //粒界の易動度
double fij[N][N]; //粒界移動の駆動力
int phinum;

int i, j, k, l, ii, jj, kk, ll, it; //整数
int ip, im, jp, jm, kp, km;         //整数
int n1, n2, n3;                     //整数

int istep = 0;
// int n000;		//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
int nstep;               //計算カウント数の最大値（計算終了カウント）
double dtime, L, dx;     // L計算領域の一辺の長さ(nm), 差分プロック１辺の長さ(m)
double M0;               //粒界の易動度
double W0;               //ペナルティー項の係数
double A0;               //勾配エネルギー係数
double F0;               //粒界移動の駆動力
double temp;             //温度
double sum1, sum2, sum3; //各種の和の作業変数
double pddtt;            //フェーズフィールドの時間変化率

double gamma0; //粒界エネルギ密度
double delta;  //粒界幅（差分ブロック数にて表現）
double mobi;   //粒界の易動度
double vm0;    //モル体積

int x11, y11, z11, x1h[10], y1h[10], z1h[10]; //初期核の座標
double t, r0, r;

//******* メインプログラム ******************************************
int main(int argc, char *argv[])
{
    nstep = 100;
    dtime = 5.0;
    temp = 1000.0;
    L = 2000.0;
    vm0 = 7.0e-6;
    delta = 7.0;
    mobi = 1.0;

    dx = L / 100 * 1.0e-9;               //差分プロック１辺の長さ(m)
    gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
    A0 = 8.0 * delta * gamma0 / PI / PI; //勾配エネルギー係数[式(4.40)]
    W0 = 4.0 * gamma0 / delta;           //ペナルティー項の係数[式(4.40)]
    M0 = mobi * PI * PI / (8.0 * delta); //粒界の易動度[式(4.40)]
    F0 = 80.0 / RR / temp;               //粒界移動の駆動力

    for (ii = 1; ii <= nm; ii++)
    {
        for (jj = 1; jj <= nm; jj++)
        {
            wij[ii][jj] = W0;
            aij[ii][jj] = A0;
            mij[ii][jj] = M0;
            fij[ii][jj] = 0.0;
            if ((ii == nm) || (jj == nm))
            {
                fij[ii][jj] = F0;
            }
            if (ii > jj)
            {
                fij[ii][jj] = -fij[ii][jj];
            }
            if (ii == jj)
            {
                wij[ii][jj] = 0.0;
                aij[ii][jj] = 0.0;
                mij[ii][jj] = 0.0;
                fij[ii][jj] = 0.0;
            }
        }
    }

    int taskid,
        numworkers,
        numtasks,
        rows, offset,
        dest, source,
        up, down,
        msgtype,
        rc, start, end,
        ix, iy, iz, it;

    MPI_Status status;

    // Allocate taskid to each core (cores = tasks = master + workers)
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    numworkers = numtasks - 1;
    rows = NDX / numworkers;

    /************************* master code *******************************/
    if (taskid == MASTER)
    {
        clock_t start_t, end_t, total_t;
        start_t = clock();
        if (NDX % numworkers != 0)
        {
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(1);
        }

        double(*phi)[N][NDX][NDY][NDZ] = malloc(sizeof(*phi));
        double(*intphi)[NDX][NDY][NDZ] = malloc(sizeof(*intphi));

        for (i = 0; i <= ndmx; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    for (ii = 1; ii <= nm - 1; ii++)
                    {
                        (*phi)[ii][i][j][k] = 0.0;
                    }
                    (*phi)[nm][i][j][k] = 1.0; // nm番目のフェーズフィールドを１に初期化
                }
            }
        }

        r0 = 10.0;
        for (ii = 1; ii <= nm - 1; ii++)
        {
            // x1=x1h[ii]; y1=y1h[ii];
            x11 = rand() % NDX;
            y11 = rand() % NDY; //初期核の位置
            z11 = rand() % NDZ;

            for (i = 0; i <= ndmx; i++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    for (k = 0; k <= ndmy; k++)
                    {
                        r = sqrt(((i - x11)) * (i - x11) + (j - y11) * (j - y11) + (k - z11) * (k - z11));
                        if (r <= r0)
                        {
                            (*phi)[ii][i][j][k] = 1.0;
                            (*phi)[nm][i][j][k] = 0.0;
                        } //初期核位置のフェーズフィールドを設定
                    }
                }
            }
        }

        offset = 0;
        // Send to workers
        for (i = 1; i <= numworkers; i++)
        {
            dest = i;
            if (dest == 1)
                up = NONE;
            else
                up = dest - 1;
            if (dest == numworkers)
                down = NONE;
            else
                down = dest + 1;

            MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&up, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&down, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            //// send phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Send(&(*phi)[ii][offset], rows * NDY * NDZ, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
            }

            offset = offset + rows;
        }

        // Receive from workers
        for (i = 1; i <= numworkers; i++)
        {
            source = i;
            msgtype = DONE;
            MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
            //// receive phase fields
            MPI_Recv(&(*intphi)[offset], rows * NDY * NDZ, MPI_DOUBLE, source,
                     msgtype, MPI_COMM_WORLD, &status);
        }

        // FILE *stream;
        // char buffer[30];
        // sprintf(buffer, "data/2d%d.csv", nstep);
        // stream = fopen(buffer, "a");

        // for (int i = 0; i <= ndmx; i++)
        // {
        //     for (int j = 0; j <= ndmy; j++)
        //     {
        //         fprintf(stream, "%e\n", intphi[i][j]);
        //     }
        // }

        // fclose(stream);

        FILE *stream;
        char buffer[30];
        sprintf(buffer, "data/2d%d.vtk", nstep);
        stream = fopen(buffer, "a");

        fprintf(stream, "# vtk DataFile Version 1.0\n");
        fprintf(stream, "phi_%d.vtk\n", nstep);
        fprintf(stream, "ASCII\n");
        fprintf(stream, "DATASET STRUCTURED_POINTS\n");
        fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
        fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
        fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
        fprintf(stream, "\n");
        fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
        fprintf(stream, "SCALARS scalars float\n");
        fprintf(stream, "LOOKUP_TABLE default\n");

        for (k = 0; k <= ndmz; k++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (i = 0; i <= ndmx; i++)
                {
                    // fprintf(streamc0, "%e\n", phi[1][i][j][k]);
                    fprintf(stream, "%e\n", (*intphi)[i][j][k]);
                }
            }
        }
        fclose(stream);

        end_t = clock();
        total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
        printf("Total time taken: %lu secs\n", total_t);
        printf("Exiting of the program...\n");

        MPI_Finalize();
    }

    /************************* workers code **********************************/
    if (taskid != MASTER)
    {
        double(*phi)[N][rows + 2][NDY][NDZ] = malloc(sizeof(*phi));
        double(*phi2)[N][rows + 2][NDY][NDZ] = malloc(sizeof(*phi2));
        int(*phiNum)[rows + 2][NDY][NDZ] = malloc(sizeof(*phiNum));
        int(*phiIdx)[N + 1][rows + 2][NDY][NDZ] = malloc(sizeof(*phiIdx));
        double(*intphi)[rows + 2][NDY][NDZ] = malloc(sizeof(*intphi));

        // Receive from master
        source = MASTER;
        msgtype = BEGIN;
        MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        //// receive phase fields
        for (ii = 1; ii <= nm; ii++)
        {
            MPI_Recv(&(*phi)[ii][1], rows * NDY * NDZ, MPI_DOUBLE, source, msgtype, MPI_COMM_WORLD, &status);
        }

    start:;

        // Communicate with neighor workers before computation
        if (up != NONE)
        {
            //// send up boundaries of phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Send(&(*phi)[ii][1][0], NDY * NDZ / 2, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
                MPI_Send(&(*phi)[ii][1][NDY / 2], NDY * NDZ / 2, MPI_DOUBLE, up, DTAG, MPI_COMM_WORLD);
            }

            source = up;
            msgtype = UTAG;
            //// receive up boundaries of phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Recv(&(*phi)[ii][0][0], NDY * NDZ / 2, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&(*phi)[ii][0][NDY / 2], NDY * NDZ / 2, MPI_DOUBLE, source,
                         msgtype, MPI_COMM_WORLD, &status);
            }
        }
        if (down != NONE)
        {
            //// send down boundaries of phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Send(&(*phi)[ii][rows][0], NDY * NDZ / 2, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
                MPI_Send(&(*phi)[ii][rows][NDY / 2], NDY * NDZ / 2, MPI_DOUBLE, down,
                         UTAG, MPI_COMM_WORLD);
            }

            source = down;
            msgtype = DTAG;
            //// receive down boundaries of phase fields
            for (ii = 1; ii <= nm; ii++)
            {
                MPI_Recv(&(*phi)[ii][rows + 1][0], NDY * NDZ / 2, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(&(*phi)[ii][rows + 1][NDY / 2], NDY * NDZ / 2, MPI_DOUBLE, source, msgtype,
                         MPI_COMM_WORLD, &status);
            }
        }

        // Compute after sending and receiving data
        start = 1;
        end = rows;

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    kp = k + 1;
                    km = k - 1;
                    if (up == NONE && i == 1)
                    {
                        im = 1;
                    }
                    if (down == NONE && i == rows)
                    {
                        ip = rows;
                    }
                    if (j == ndmy)
                    {
                        jp = 0;
                    }
                    if (j == 0)
                    {
                        jm = ndmy;
                    }
                    if (k == ndmz)
                    {
                        kp = 0;
                    }
                    if (k == 0)
                    {
                        km = ndmz;
                    }

                    //--- 位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数---
                    phinum = 0;
                    for (ii = 1; ii <= nm; ii++)
                    {
                        if (((*phi)[ii][i][j][k] > 0.0) ||
                            (((*phi)[ii][i][j][k] == 0.0) && ((*phi)[ii][ip][j][k] > 0.0) ||
                             ((*phi)[ii][im][j][k] > 0.0) ||
                             ((*phi)[ii][i][jp][k] > 0.0) ||
                             ((*phi)[ii][i][jm][k] > 0.0) ||
                             ((*phi)[ii][i][j][kp] > 0.0) ||
                             ((*phi)[ii][i][j][km] > 0.0)))
                        {
                            phinum++;
                            (*phiIdx)[phinum][i][j][k] = ii;
                            // printf("%d  ", n00);
                        }
                    }
                    (*phiNum)[i][j][k] = phinum;
                }
            }
        }

        // Evolution Equations
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    kp = k + 1;
                    km = k - 1;
                    if (up == NONE && i == 1)
                    {
                        im = 1;
                    }
                    if (down == NONE && i == rows)
                    {
                        ip = rows;
                    }
                    if (j == ndmy)
                    {
                        jp = 0;
                    }
                    if (j == 0)
                    {
                        jm = ndmy;
                    }
                    if (k == ndmz)
                    {
                        kp = 0;
                    }
                    if (k == 0)
                    {
                        km = ndmz;
                    }

                    for (n1 = 1; n1 <= (*phiNum)[i][j][k]; n1++)
                    {
                        ii = (*phiIdx)[n1][i][j][k];
                        pddtt = 0.0;
                        for (n2 = 1; n2 <= (*phiNum)[i][j][k]; n2++)
                        {
                            jj = (*phiIdx)[n2][i][j][k];
                            sum1 = 0.0;
                            for (n3 = 1; n3 <= (*phiNum)[i][j][k]; n3++)
                            {
                                kk = (*phiIdx)[n3][i][j][k];
                                sum1 += 0.5 * (aij[ii][kk] - aij[jj][kk]) * ((*phi)[kk][ip][j][k] + (*phi)[kk][im][j][k] + (*phi)[kk][i][jp][k] + (*phi)[kk][i][jm][k] + (*phi)[kk][i][j][kp] + (*phi)[kk][i][j][km] - 6.0 * (*phi)[kk][i][j][k]) + (wij[ii][kk] - wij[jj][kk]) * (*phi)[kk][i][j][k]; //[式(4.31)の一部]
                            }
                            pddtt += -2.0 * mij[ii][jj] / (double)((*phiNum)[i][j][k]) * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt((*phi)[ii][i][j][k] * (*phi)[jj][i][j][k]));
                            //フェーズフィールドの発展方程式[式(4.31)]
                        }
                        (*phi2)[ii][i][j][k] = (*phi)[ii][i][j][k] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
                        if ((*phi2)[ii][i][j][k] >= 1.0)
                        {
                            (*phi2)[ii][i][j][k] = 1.0;
                        } //フェーズフィールドの変域補正
                        if ((*phi2)[ii][i][j][k] <= 0.0)
                        {
                            (*phi2)[ii][i][j][k] = 0.0;
                        }
                    }
                } // k
            }     // j
        }         // i

        for (ii = 1; ii <= nm; ii++)
        {
            for (i = start; i <= end; i++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    for (k = 0; k <= ndmz; k++)
                    {
                        (*phi)[ii][i][j][k] = (*phi2)[ii][i][j][k];
                    }
                }
            }
        }

        //
        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    sum1 = 0.0;
                    for (ii = 1; ii <= nm; ii++)
                    {
                        sum1 += (*phi)[ii][i][j][k];
                    }
                    for (ii = 1; ii <= nm; ii++)
                    {
                        (*phi)[ii][i][j][k] = (*phi)[ii][i][j][k] / sum1;
                    }
                }
            }
        }

        istep = istep + 1.;
        if (istep < nstep)
        {
            goto start;
        }
    end:;

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    for (ii = 1; ii <= nm; ii++)
                    {
                        (*intphi)[i][j][k] += (*phi)[ii][i][j][k] * (*phi)[ii][i][j][k];
                    }
                }
            }
        }
        // Send final result to master
        MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        //// send phase fields
        MPI_Send(&(*intphi)[1], rows * NDY * NDZ, MPI_DOUBLE, MASTER, DONE,
                 MPI_COMM_WORLD);
        MPI_Finalize();
    }
    return 0;
}
