// LAB1.cpp: определяет точку входа для консольного приложения. 
// 

//#include "stdafx.h" 
#include <stdio.h> 
#include <math.h> 
#include <stdlib.h> 
#include <ctime> 
#include "mpi.h" 
#include <iostream>
#include <windows.h>

enum ConsoleColor
{
	Black = 0,
	Blue = 1,
	Green = 2,
	Cyan = 3,
	Red = 4,
	Magenta = 5,
	Brown = 6,
	LightGray = 7,
	DarkGray = 8,
	LightBlue = 9,
	LightGreen = 10,
	LightCyan = 11,
	LightRed = 12,
	LightMagenta = 13,
	Yellow = 14,
	White = 15
};

void SetColor(ConsoleColor text, ConsoleColor background)
{
	HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hStdOut, (WORD)((background << 4) | text));
}


void Laplas(float *matr, int N, int &_iter)
{
	float eps = 0.1;

	float dmax, temp, dm;
	int elem_count = (N + 1)*(N + 1);
	float *u = new float[elem_count];
	for (int i = 0; i < elem_count; i++)
		u[i] = matr[i];

	float ctrl;
	do {
		dmax = 0; // максимальное изменение значений u
		for (int i = 1; i < (N); i++)
			for (int j = 1; j < (N); j++) {

				temp = u[i * (N) + j];
				ctrl = 0.25*(u[(i - 1)*(N)+j] + u[(i + 1)*(N)+j] + u[i*(N)+(j - 1)] + u[i*(N)+(j + 1)]);
				u[i*(N)+j] = 0.25*(u[(i - 1)*(N)+j] + u[(i + 1)*(N)+j] + u[i*(N)+(j - 1)] + u[i*(N)+(j + 1)]);
				dm = fabs(temp - u[i*(N)+j]); // abs
				if (dmax < dm) dmax = dm;
			}
		_iter++;
	} while (dmax > eps);
}

int main(int argc, char* argv[])
{

	int proc_num, proc_rank;

	int chunk_size;
	double parallel_start, parallel_end;
	double linear_start, linear_end;
	int iter = 0, yo = 0;

	float eps = 0.1;

	int N = 499;

	float ftemp, dm;
	int elem_count = (N + 1)*(N + 1);

	float *u = new float[elem_count];
	float *peredacha = 0;

	//---------------------------------------Начало---------------------------------------
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
	float dmax = 0;

	int NB = sqrt(proc_num); //!!!!!!!!!!!!!!!!!!!!!!

	chunk_size = (N + 1) / NB;

	// Создание bufs для граничных точек
	float *top = new float[chunk_size], *left = new float[chunk_size], *bot = new float[chunk_size], *right = new float[chunk_size];
	
	if (proc_rank == 0)
	{
		std::cout << "Chunk_Size = " << chunk_size << std::endl;
		peredacha = new float[elem_count];
		for (int i = 0; i < elem_count; i++)
			u[i] = rand() % (N * 2) - N;

		srand(MPI_Wtime() * 1000);

		for (int i = 0; i < elem_count; i++)
			peredacha[i] = u[i];

		parallel_start = MPI_Wtime();

	}

	float *chunk = new float[chunk_size*chunk_size];

	MPI_Datatype block; // Особый тип данных
	MPI_Type_vector(
		chunk_size, // Число блоков
		chunk_size, // Число элементов базового типа в каждом блоке
		N + 1,      // Шаг между началами соседних блоков, измеренный числом элементов базового типа
		MPI_FLOAT,  // Базовый тип данных
		&block);    // Новый производный тип данных
	MPI_Type_commit(&block); // Регисрация нового произвольного типа

	MPI_Status status;

	if (proc_rank == 0) // Заносим в chunk массив блоков на 0 процесс
		// posilaem i prinimaem chast v 0
		MPI_Sendrecv(
		&u[0], 1, block, 0, 8, // Параметры передаваемого сообщения
		chunk, chunk_size*chunk_size, MPI_FLOAT, 0, 8, // Параметры принимаемого сообщения
		MPI_COMM_WORLD, &status); 

	if (proc_rank == 0)	// Рассылаем по блокам на процессы																								// rassilaem vsem ostalnim
		for (int s = 1; s < proc_num; s++)
			MPI_Send(&u[((s / NB)*(N + 1) + (s % NB))*chunk_size], 1, block, s, 9, MPI_COMM_WORLD);

	if (proc_rank != 0) // Принимаем
		MPI_Recv(chunk, chunk_size*chunk_size, MPI_FLOAT, 0, 9, MPI_COMM_WORLD, &status);


	float i1, i2, j1, j2;
	int ek = 0, el = 0, lk = 0, ll = 0;

	if (proc_rank / NB == 0) ek = 1;
	if (proc_rank % NB == 0) el = 1;
	if (proc_rank / NB == NB - 1) lk = -1;
	if (proc_rank % NB == NB - 1) ll = -1;
	float fer;

	do {
		dmax = 0;

		// Получение граничных узлов

		if (proc_rank / NB != 0) { // Строка не нулевая

			// Получение данных от верхнего процессора
			MPI_Status status;

			MPI_Recv(top, chunk_size, MPI_FLOAT, proc_rank - NB, 10, MPI_COMM_WORLD, &status);

			// Пересылка данных verhnemy процессору
			// Verhnyaa строка
			float *temp = new float[chunk_size];
			for (int i = 0; i < chunk_size; i++)
				temp[i] = chunk[i]; // Посылка строки
			MPI_Send(temp, chunk_size, MPI_FLOAT, proc_rank - NB, 1, MPI_COMM_WORLD);

		}

		if (proc_rank % NB != 0) { // столбец не нулевой
			//prinyat ot levogo
			MPI_Status status;
			int k = MPI_Recv(left, chunk_size, MPI_FLOAT, proc_rank - 1, 2, MPI_COMM_WORLD, &status);
			// пересылка pravogo данных levomy процессору
			float *temp = new float[chunk_size];
			for (int i = 0; i < chunk_size; i++)
				temp[i] = chunk[i * chunk_size]; // Посылка столбца
			MPI_Send(temp, chunk_size, MPI_FLOAT, proc_rank - 1, 3, MPI_COMM_WORLD);
		}

		// пересылка граничных узлов
		if (proc_rank / NB != NB - 1) {
			// пересылка данных verhnemu процессору
			float *temp = new float[chunk_size];
			for (int i = 0; i < chunk_size; i++)
				temp[i] = chunk[chunk_size*(chunk_size - 1) + i];
			MPI_Send(temp, chunk_size, MPI_FLOAT, proc_rank + NB, 10, MPI_COMM_WORLD);

			MPI_Status status;
			int k = MPI_Recv(bot, chunk_size, MPI_FLOAT, proc_rank + NB, 1, MPI_COMM_WORLD, &status);
		}

		if (proc_rank % NB != NB - 1) {
			//pravomu pravoe
			float *temp = new float[chunk_size];
			for (int i = 0; i<chunk_size; i++)
				temp[i] = chunk[i*chunk_size + chunk_size - 1];
			MPI_Send(temp, chunk_size, MPI_FLOAT, proc_rank + 1, 2, MPI_COMM_WORLD);

			MPI_Status status;
			int k = MPI_Recv(right, chunk_size, MPI_FLOAT, proc_rank + 1, 3, MPI_COMM_WORLD, &status);
		}


		// Обработка блока с оценкой погрешности dmax

		for (int i = ek; i < chunk_size + lk; i++)
			for (int j = el; j < chunk_size + ll; j++) {
				ftemp = chunk[i * (chunk_size) + j];

				// Вычисление соседей
				if (i - 1 >= 0) i1 = chunk[(i - 1) * (chunk_size) + j];
				else i1 = top[j];
				if (i + 1 < chunk_size) i2 = chunk[(i + 1) * (chunk_size) + j];
				else i2 = bot[j];
				if (j - 1 >= 0) j2 = chunk[i * (chunk_size) + (j - 1)];
				else j2 = left[i];
				if (j + 1 < chunk_size) j1 = chunk[i * (chunk_size) + (j + 1)];
				else j1 = right[i];

				chunk[i * (chunk_size) + j] = 0.25 * (i1 + i2 + j1 + j2);

				dm = fabs(ftemp - chunk[i*(chunk_size)+j]); // abs
				if (dmax < dm) dmax = dm;
			}
		iter++;
		// синхронизация и рассылка погрешности dmax

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Reduce(&dmax, &fer, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		dmax = fer;
		MPI_Bcast(&dmax, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	} while (dmax > eps);

	if (proc_rank == 0)
	{
		parallel_end = MPI_Wtime();

		std::cout << std::endl << "-----Parallel version-----" << std::endl;
		std::cout << "Time: " << parallel_end - parallel_start << " sec" << std::endl;
		std::cout << "The number of operation: " << iter << std::endl << std::endl;


		linear_start = MPI_Wtime();

		Laplas(peredacha, N, iter);

		linear_end = MPI_Wtime();

		std::cout << std::endl << "-----Linear version-----" << std::endl;
		std::cout << "Time: " << linear_end - linear_start << " sec" << std::endl;
		std::cout << "The number of operation: " << iter << std::endl << std::endl << std::endl;

		std::cout << "The time difference: "
			<< linear_end - linear_start - parallel_end + parallel_start
			<< std::endl;
		std::cout << "Acceleration: "
			<< (linear_end - linear_start) / (parallel_end - parallel_start)
			<< std::endl << std::endl;

		int flag = 1;
		for (int i = 0; i < chunk_size*chunk_size; i++)
			if (chunk[i] != peredacha[i]) flag = 0;

		std::cout << "Check for correctness: ";
		if (flag) { std::cout << "[";  SetColor(Red, Black); std::cout << " Faild "; SetColor(White, Black); std::cout << "]" << std::endl << std::endl; }
		else { std::cout << "[";  SetColor(LightGreen, Black); std::cout << " Correctly "; SetColor(White, Black); std::cout << "]" << std::endl << std::endl; }
	}

	MPI_Finalize();

	return 0;
}