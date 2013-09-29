#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mpi.h"

// Problem parameters.
int screen_width, screen_height;
int wire_width, wire_height;
int wire_shift_x, wire_shift_y;
int iteration_limit;
double wire_voltage;

// Algorithm parameters.
double v_top, v_right, v_bottom, v_left, v_x;

// Data rows.
// Each segment border shares with other adjacent segments.
double* data;

// Process data size;
int segment_size_x, segment_size_y;

// Board size.
int process_board_x, process_board_y;

// Position on board for current process_number.
int process_x, process_y;

// MPI specific definitions.
int available_processes, process_number;
MPI_Datatype MPI_column;

/*
#ifdef DEBUG
#define d(format, args...)
#else
#define d(format, args...) _debug
#endif
*/

// former _debug
void d(const char*, ...);

void print_board();

double initial_voltage(const int, const int);
int is_wire(const int, const int);
int is_outside(const int, const int);

double calculation(const int);
double* get_point_address(int, int);

int top_point(const int, const int);
int right_point(const int, const int);
int bottom_point(const int, const int);
int left_point(const int, const int);

int top_segment();
int right_segment();
int bottom_segment();
int left_segment();

int file_exists(const char*);

void finalize();
void communicate(int);

int main(int argc, char* argv[]) {
	int x, y;
	char buffer[255];
	FILE* file;
	int iteration = 0;
	int progress = -1;
	int add_headers = 0;
	double timer;

	if (argc < 7) {
		printf("Usage: %s screen_width screen_height wire_width wire_height wire_mili_voltage iteration_limit\n", argv[0]);
		return 1;
	}

	// Screen size.
	screen_width = atoi(argv[1]);
	screen_height = atoi(argv[2]);

	// Wire inside screen size.
	wire_width = atoi(argv[3]);
	wire_height = atoi(argv[4]);

	if (screen_width <= wire_width) {
		screen_width = wire_width;
		wire_width = atoi(argv[1]);
	}

	if (screen_height <= wire_height) {
		screen_height = wire_height;
		wire_height = atoi(argv[2]);
	}

	// Wire constant voltage; input as mV thus "/ 1000" part.
	wire_voltage = atoi(argv[5]) / 1000;

	// Maximum number of iterations.
	iteration_limit = atoi(argv[6]);

	d("Initiate MPI world%s", ".");
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &available_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_number);

	if (available_processes < 2) {
		printf("%s can be used only with at least two nodes\n", argv[0]);
		finalize();
		return 0;
	}

	wire_shift_x = (screen_width - wire_width) / 2;
	wire_shift_y = (screen_height - wire_height) / 2;

	process_board_x = sqrt(available_processes);
	process_board_y = available_processes / process_board_x;
	process_x = process_number % process_board_x;
	process_y = process_number / process_board_x;

	if (0 == process_number) {
		print_board();
	}

	printf("My rank: %d.\nMy location: %d:%d (board size: %d:%d)\n", process_number, process_x, process_y, process_board_x, process_board_y);
	if (process_number >= process_board_x * process_board_y) {
		printf("I'm out of board. Destroy MPI world.\n");
		finalize();
		return 0;
	}

	d("Initiate data storage%s", ".");
	segment_size_x = screen_width / process_board_x + 2;
	segment_size_y = screen_height / process_board_y + 2;
	data = (double*) malloc(sizeof(double) * segment_size_x * segment_size_y);

	for (y = 0; y < segment_size_y; ++y) {
		for (x = 0; x < segment_size_x; ++x) {
			data[y * segment_size_x + x] = initial_voltage(x, y);
		}
	}

	d("Register custom MPI type%s", ".");
	MPI_Type_vector(segment_size_y - 2, 1, segment_size_x, MPI_DOUBLE, &MPI_column);
	MPI_Type_commit(&MPI_column);

	timer = MPI_Wtime();
	communicate(iteration);
	while (iteration < iteration_limit) {
		d("Iterations start%s", ".");
		calculation(iteration);
		communicate(iteration);
		++iteration;
		if (process_number == 0 && 20 * iteration / iteration_limit > progress) {
			progress = 20 * iteration / iteration_limit;
			printf("Progress: %d%% (%d of %d)\n", 5 * progress, iteration, iteration_limit);
		}
	}
	timer = MPI_Wtime() - timer;

	// Dump data from each process.
	sprintf(buffer, "data-%dx%d.txt", process_x, process_y);
	file = fopen(buffer, "w");
	for (y = 1; y < segment_size_y - 1; ++y) {
		for (x = 1; x < segment_size_x - 1; ++x) {
			fprintf(file, "%f\t%f\t%f\n",
					1.0 * (process_x * (segment_size_x - 2) + x) / screen_width,
					1.0 * (process_y * (segment_size_y - 2) + y) / screen_height,
					data[y * segment_size_x + x]
			);
		}
		fprintf(file, "\n");
	}
	fclose(file);

	// Log program response.
	if (process_number == 0) {
		printf("Board size: %dx%d, iterations limit: %d, used processes: %d, elapsed time: %f\n", screen_width, screen_height, iteration_limit, available_processes, timer);

		if (file_exists("results.txt") == 0) {
			add_headers = 1;
		}

		file = fopen("results.txt", "a+");
		if (add_headers) {
			fprintf(file, "W\tH\tmax\tproc\tT\n");
		}
		fprintf(file, "%d\t%d\t%d\t%d\t%f\n", screen_width, screen_height, iteration_limit, available_processes, timer);
		fclose(file);
	}

	finalize();

	return 0;
}

// former _debug
void d(const char* format, ...) {
	va_list args;
	va_start(args, format);

/**/
	char _format[500];
	sprintf(_format, "%s:%d (%s): %s", __FILE__, __LINE__, __FUNCTION__, format);
	sprintf(_format, "%d (%d:%d): %s\n", process_number, process_x, process_y, _format);
	vprintf(_format, args);
/**/
	va_end(args);
}

void print_board() {
	printf("Board layout:\n");
	for (int i = 0; i < process_board_x; ++i) {
		if (0 == i) {
			printf("-");
			for (int j = 0; j < process_board_y; ++j) {
				printf("-----");
			}
			printf("\n");
		}
		printf("|");
		for (int j = 0; j < process_board_y; ++j) {
			printf(" %02d |", j + i * process_board_y);
		}
		printf("\n-");
		for (int j = 0; j < process_board_y; ++j) {
			printf("-----");
		}
		printf("\n");
	}
}

double initial_voltage(const int x, const int y) {
	if (is_wire(x, y)) {
		return wire_voltage;
	}

	return 0.0;
}

int is_wire(const int x, const int y) {
	int board_x, board_y;

	board_x = process_x * (segment_size_x - 2) + x - 1;
	board_y = process_y * (segment_size_y - 2) + y - 1;

	if (board_x < wire_shift_x
		|| board_y < wire_shift_y
		|| board_x >= wire_shift_x + wire_width
		|| board_y >= wire_shift_y + wire_height
	) {
		return 0;
	}

	return 1;
}

int is_outside(const int x, const int y) {
	int board_x, board_y;

	board_x = process_x * (segment_size_x - 2) + x - 1;
	board_y = process_y * (segment_size_y - 2) + y - 1;

	if (board_x < 0
		|| board_y < 0
		|| board_x >= screen_width
		|| board_y >= screen_height
	) {
		return 1;
	}

	return 0;
}

// Calculate one step and max corrections
double calculation(const int iterate) {
	// Vt + Vr + Vb + Vl - 4Vx = 0
	int y, x;
	double maxStepSize, pde, stepSize;
	d("Calculation step (iterate: %d)", iterate);
	maxStepSize = 0;

	for (y = 1; y < segment_size_y - 1; ++y) {
		for (x = 1; x < segment_size_x - 1; ++x) {
			// We're using checkerboard to calculate values.
			if ((y + x + iterate) % 2 == 0
				|| is_outside(x, y)
				|| is_wire(x, y)
			) {
				continue;
			}

			pde = data[bottom_point(x, y)]
				+ data[top_point(x, y)]
				+ data[left_point(x, y)]
				+ data[right_point(x, y)];
			pde *= .25;

			stepSize = abs(pde - data[y * segment_size_x + x]);
			if (stepSize > maxStepSize) {
				maxStepSize = stepSize;
			}

			data[y * segment_size_x + x] = pde;
		}
	}

	return maxStepSize;
}

// Get point data; wrap coordinates to fit segment box.
double* get_point_address(int x, int y) {
	if (x < 0) {
		x += segment_size_x;
	} else {
		while (x >= segment_size_x) {
			x -= segment_size_x;
		}
	}
	if (y < 0) {
		y += segment_size_y;
	} else {
		while (y >= segment_size_y) {
			y -= segment_size_y;
		}
	}

	return data + y * segment_size_x + x;
}

int left_point(const int x, const int y) {
	if (x > 0) {
		return y * segment_size_x + x - 1;
	}

	return -1;
}

int right_point(const int x, const int y) {
	if (1 + x < segment_size_x) {
		return y * segment_size_x + x + 1;
	}

	return -1;
}

int top_point(const int x, const int y) {
	if (y > 0) {
		return (y - 1) * segment_size_x + x;
	}

	return -1;
}

int bottom_point(const int x, const int y) {
	if (1 + y < segment_size_y) {
		return (y + 1) * segment_size_x + x;
	}

	return -1;
}

int normalize_segment(int point) {
	return (-1 == point) ? point : ceil(1.0 * point * (process_board_x * process_board_y) / (screen_width * screen_height));
}

int left_segment() {
  int point = left_point(process_x, process_y);
	int segment = normalize_segment(point);

  return segment;
}

int right_segment() {
  int point = right_point(process_x, process_y);
	int segment = normalize_segment(point);

  return segment;

	return (-1 == point) ? point : ceil(1.0 * point / (process_board_x * process_board_y));
}

int top_segment() {
  int point = top_point(process_x, process_y);
	int segment = normalize_segment(point);

  return segment;

	return (-1 == point) ? point : ceil(1.0 * point / (process_board_x * process_board_y));
}

int bottom_segment() {
  int point = bottom_point(process_x, process_y);
	int segment = normalize_segment(point);

  return segment;

	return (-1 == point) ? point : ceil(1.0 * point / (process_board_x * process_board_y));
}

int file_exists(const char* filename) {
	FILE* file = fopen(filename, "r");
	if (file) {
		fclose(file);
		return 1;
	}
	return 0;
}

// Communicate with adjacent segments.
void communicate(int iteration) {
	MPI_Request left_request, right_request, top_request, bottom_request;
	MPI_Status status;

	if (top_segment() > -1) {
		d("Sending top row to: %d", top_segment());
		MPI_Isend(get_point_address(1, 1), segment_size_x - 2, MPI_DOUBLE, top_segment(), iteration, MPI_COMM_WORLD, &top_request);
	}
	if (right_segment() > -1) {
		d("Sending right column to: %d", right_segment());
		MPI_Isend(get_point_address(-2, 1), 1, MPI_column, right_segment(), iteration, MPI_COMM_WORLD, &right_request);
	}
	if (bottom_segment() > -1) {
		d("Sending bottom row to: %d", bottom_segment());
		MPI_Isend(get_point_address(1, -2), segment_size_x - 2, MPI_DOUBLE, bottom_segment(), iteration, MPI_COMM_WORLD, &bottom_request);
	}
	if (left_segment() > -1) {
		d("Sending left column to: %d", left_segment());
		MPI_Isend(get_point_address(1, 1), 1, MPI_column, left_segment(), iteration, MPI_COMM_WORLD, &left_request);
	}


	if (top_segment() > -1) {
		d("Receiving top row from: %d", top_segment());
		MPI_Recv(get_point_address(1, 0), segment_size_x - 2, MPI_DOUBLE, top_segment(), iteration, MPI_COMM_WORLD, &status);
	}
	if (right_segment() > -1) {
		d("Receiving right column from: %d", right_segment());
		MPI_Recv(get_point_address(-1, 1), 1, MPI_column, right_segment(), iteration, MPI_COMM_WORLD, &status);
	}
	if (bottom_segment() > -1) {
		d("Receiving bottom row from: %d", bottom_segment());
		MPI_Recv(get_point_address(1, -1), segment_size_x - 2, MPI_DOUBLE, bottom_segment(), iteration, MPI_COMM_WORLD, &status);
	}
	if (left_segment() > -1) {
		d("Receiving left column from: %d", left_segment());
		MPI_Recv(get_point_address(0, 1), 1, MPI_column, left_segment(), iteration, MPI_COMM_WORLD, &status);
	}

	if (top_segment() > -1) {
		d("Wait for top MPI requests (%d).", top_segment());
		MPI_Wait(&top_request, &status);
	}
	if (right_segment() > -1) {
		d("Wait for right MPI requests (%d).", right_segment());
		MPI_Wait(&right_request, &status);
	}
	if (bottom_segment() > -1) {
		d("Wait for bottom MPI requests (%d).", bottom_segment());
		MPI_Wait(&bottom_request, &status);
	}
	if (left_segment() > -1) {
		d("Wait for left MPI requests (%d).", left_segment());
		MPI_Wait(&left_request, &status);
	}
}

void finalize() {
	MPI_Finalize();
}
