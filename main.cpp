#include "inmost.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace INMOST;

namespace problem
{
    double u_exact(double x, double y)
    {
        return std::sin(4.0 * x) * std::cos(3.0 * y);
    }

    double rhs(double x, double y)
    {
        return 25.0 * u_exact(x, y);
    }

    double boundary_bottom(double x)
    {
        return u_exact(x, 0.0);
    }

    double boundary_top(double x)
    {
        return u_exact(x, 1.0);
    }

    double boundary_left(double y)
    {
        return u_exact(0.0, y);
    }

    double boundary_right(double y)
    {
        return u_exact(1.0, y);
    }
}

unsigned node_id(unsigned row, unsigned col, unsigned grid_size)
{
    return row * grid_size + col;
}

void assemble_matrix(Sparse::Matrix &matrix, unsigned grid_size)
{
    const unsigned total_unknowns = grid_size * grid_size;

    matrix.SetInterval(0, total_unknowns);

    for (unsigned row = 0; row < grid_size; ++row)
    {
        for (unsigned col = 0; col < grid_size; ++col)
        {
            unsigned id = node_id(row, col, grid_size);

            matrix[id][id] = 4.0;

            if (col > 0)
                matrix[id][node_id(row, col - 1, grid_size)] = -1.0;

            if (col + 1 < grid_size)
                matrix[id][node_id(row, col + 1, grid_size)] = -1.0;

            if (row > 0)
                matrix[id][node_id(row - 1, col, grid_size)] = -1.0;

            if (row + 1 < grid_size)
                matrix[id][node_id(row + 1, col, grid_size)] = -1.0;
        }
    }
}

void assemble_rhs(Sparse::Vector &rhs_vec, unsigned grid_size, double step)
{
    const unsigned total_unknowns = grid_size * grid_size;

    rhs_vec.SetInterval(0, total_unknowns);

    for (unsigned row = 0; row < grid_size; ++row)
    {
        double y = (row + 1) * step;

        for (unsigned col = 0; col < grid_size; ++col)
        {
            double x = (col + 1) * step;
            unsigned id = node_id(row, col, grid_size);

            rhs_vec[id] = problem::rhs(x, y) * step * step;

            if (row == 0)
                rhs_vec[id] += problem::boundary_bottom(x);

            if (row + 1 == grid_size)
                rhs_vec[id] += problem::boundary_top(x);

            if (col == 0)
                rhs_vec[id] += problem::boundary_left(y);

            if (col + 1 == grid_size)
                rhs_vec[id] += problem::boundary_right(y);
        }
    }
}

double numerical_value(
    const Sparse::Vector &solution,
    unsigned i,
    unsigned j,
    unsigned grid_size,
    double step)
{
    double x = j * step;
    double y = i * step;

    if (i == 0)
        return problem::boundary_bottom(x);

    if (i == grid_size + 1)
        return problem::boundary_top(x);

    if (j == 0)
        return problem::boundary_left(y);

    if (j == grid_size + 1)
        return problem::boundary_right(y);

    return solution[node_id(i - 1, j - 1, grid_size)];
}

void calculate_errors(
    const Sparse::Vector &solution,
    unsigned grid_size,
    double step,
    double &max_error,
    double &l2_error)
{
    max_error = 0.0;
    l2_error = 0.0;

    for (unsigned i = 0; i <= grid_size + 1; ++i)
    {
        double y = i * step;

        for (unsigned j = 0; j <= grid_size + 1; ++j)
        {
            double x = j * step;

            double exact = problem::u_exact(x, y);
            double approx = numerical_value(solution, i, j, grid_size, step);

            double difference = std::abs(approx - exact);

            max_error = std::max(max_error, difference);
            l2_error += difference * difference;
        }
    }

    l2_error = std::sqrt(l2_error * step * step);
}

int main(int argc, char **argv)
{
    Solver::Initialize(&argc, &argv);

    std::vector<unsigned> grid_sizes =
    {
        10, 20, 40, 80, 160, 320, 640, 1280
    };

    std::cout << std::setw(8)  << "h"
              << std::setw(16) << "C-norm error"
              << std::setw(16) << "L2-norm error"
              << std::setw(12) << "Iters"
              << std::setw(16) << "Iter.time"
			  << std::setw(16) << "Total.time"
              << std::endl;

    for (unsigned grid_size : grid_sizes)
    {
        const unsigned unknown_count = grid_size * grid_size;
        const double step = 1.0 / (grid_size + 1);

        Sparse::Matrix system_matrix;
        Sparse::Vector right_part;
        Sparse::Vector answer;

        answer.SetInterval(0, unknown_count);

		auto start_total = std::chrono::high_resolution_clock::now();

        assemble_matrix(system_matrix, grid_size);
        assemble_rhs(right_part, grid_size, step);

        Solver linear_solver(Solver::INNER_ILU2);

        linear_solver.SetParameter("absolute_tolerance", "1e-14");
        linear_solver.SetParameter("relative_tolerance", "1e-11");
        linear_solver.SetMatrix(system_matrix);

        bool success = linear_solver.Solve(right_part, answer);

		auto end_total = std::chrono::high_resolution_clock::now();

		double total_time = std::chrono::duration<double>(end_total - start_total).count();

        if (!success)
        {
            std::cout << "Computation failed for grid size = "
                      << grid_size << std::endl;
            continue;
        }

        double uniform_norm = 0.0;
        double integral_norm = 0.0;

        calculate_errors(
            answer,
            grid_size,
            step,
            uniform_norm,
            integral_norm
        );

        std::cout << std::setw(8) << step
                  << std::setw(16) << std::scientific << uniform_norm
                  << std::setw(16) << integral_norm
                  << std::setw(12) << linear_solver.Iterations()
                  << std::setw(16) << linear_solver.IterationsTime()
				  << std::setw(16) << total_time
                  << std::endl;
    }

    Solver::Finalize();

    return 0;
}
