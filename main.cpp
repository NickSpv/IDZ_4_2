#include "gnuplot.h"
#include "constants.h"
#include "fstream"
#include "cmath"
#include <utility>
#include <vector>

class Euler {
public:
    Euler(double _x_max, double _n, double _x_0, std::vector<double> _y_0) : x_max(_x_max), n(_n), x_0(_x_0), y_0(std::move(_y_0)) {
        h = (x_max - _x_0) / n;
    }

    void printDataInFileFor123(std::vector<std::string> file_names) {
        std::vector<std::ofstream> files(file_names.size());
        for (int i = 0; i < file_names.size(); i++) {
            files[i].open(file_names[i]);
            if (!files[i]) {
                std::cout << "Error\n";
            }
        }
        std::vector<double> y_i = y_0;
        for (int i = 0; i < y_i.size(); i++) {
            files[i] << x_0 << ' ' << y_0[i] << '\n';
            std::cout << '\n' << x_0 << ' ' << y_0[i];
        }
        for (int i = 0; i < n; i++) {
            double x_i = x_0 + i * h;
            std::vector<double> new_y_i(y_i.size());
            if (i % static_cast<int>(n / 10) == 0)
                std::cout << '\n' << x_i + h << ' ';
            for (int j = 0; j < y_i.size(); j++) {
                new_y_i[j] = y_i[j] + h * functions[j](x_i, y_i);
                files[j] << x_i + h << ' ' << new_y_i[j] << '\n';
                if(i % static_cast<int>(n / 10) == 0)
                    std::cout << new_y_i[j] << ' ';
            }
            y_i = new_y_i;
        }
        for (int i = 0; i < y_i.size(); i++) {
            files[i].close();
        }
    }

    void printDataInFileFor4(const std::string& file_name) {
        std::ofstream file;
        file.open(file_name);
        if (!file) {
            std::cout << "Error\n";
        }
        std::vector<double> y_i = y_0;
        file << x_0 << ' ' << y_i[0] << '\n';
        for (int i = 0; i < n; i++) {
            double x_i = x_0 + i * h;
            std::vector<double> new_y_i(y_i.size());
            for (int j = 0; j < y_i.size(); j++) {
                new_y_i[j] = y_i[j] + h * functions[j](x_i, y_i);
            }
            file << x_i + h << ' ' << new_y_i[0] << '\n';
            y_i = new_y_i;
        }
        file.close();
    }

    void setFunction(std::vector<double (*)(double, std::vector<double>)> _functions) {
        functions = std::move(_functions);
    }

private:
    double x_max;
    double n;
    double h;
    double x_0;
    std::vector<double> y_0;
    std::vector<double (*)(double, std::vector<double>)> functions;
};

class RungeKutta4 {
public:
    RungeKutta4(double _x_max, double _n, double _x_0, std::vector<double> _y_0) : x_max(_x_max), n(_n), x_0(_x_0), y_0(std::move(_y_0)) {
        h = (x_max - _x_0) / n;
    }

    void printDataInFileFor123(std::vector<std::string> file_names) {
        std::vector<std::ofstream> files(file_names.size());
        for (int i = 0; i < file_names.size(); i++) {
            files[i].open(file_names[i]);
            if (!files[i]) {
                std::cout << "Error\n";
            }
        }
        std::vector<double> y_i = y_0;
        for (int i = 0; i < y_i.size(); i++) {
            files[i] << x_0 << ' ' << y_0[i] << '\n';
        }
        for (int i = 0; i < n; i++) {
            double x_i = x_0 + i * h;
            std::vector<double> new_y_i(y_i.size());
            std::vector<double> K(4);
            for (int j = 0; j < y_i.size(); j++) {
                K[0] = h * functions[j](x_i, y_i);

                std::vector<double> y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 2;
                K[1] = h * functions[j](x_i + h / 2, y_i_temp);

                y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[1] / 2;
                K[2] = h * functions[j](x_i + h / 2, y_i_temp);

                y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[2];
                K[3] = h * functions[j](x_i + h, y_i_temp);

                new_y_i[j] = y_i[j] + (1 / 6.) * (K[0] + 2 * K[1] + 2 * K[2] + K[3]);
                files[j] << x_i + h << ' ' << new_y_i[j] << '\n';
            }
            y_i = new_y_i;
        }
        for (int i = 0; i < y_i.size(); i++) {
            files[i].close();
        }
    }

    void printDataInFileFor4(const std::string& file_name) {
        std::ofstream file;
        file.open(file_name);
        if (!file) {
            std::cout << "Error\n";
        }
        std::vector<double> y_i = y_0;
        file << x_0 << ' ' << y_i[0] << '\n';
        for (int i = 0; i < n; i++) {
            double x_i = x_0 + i * h;
            std::vector<double> new_y_i(y_i.size());
            std::vector<double> K(4);
            for (int j = 0; j < y_i.size(); j++) {
                K[0] = h * functions[j](x_i, y_i);

                std::vector<double> y_i_temp = y_i;
                for (auto k : y_i_temp) k += K[0] / 2;
                K[1] = h * functions[j](x_i + h / 2, y_i_temp);

                y_i_temp = y_i;
                for (auto k : y_i_temp) k += K[1] / 2;
                K[2] = h * functions[j](x_i + h / 2, y_i_temp);

                y_i_temp = y_i;
                for (auto k : y_i_temp) k += K[2];
                K[3] = h * functions[j](x_i + h, y_i_temp);

                new_y_i[j] = y_i[j] + (1 / 6.) * (K[0] + 2 * K[1] + 2 * K[2] + K[3]);
            }
            file << x_i + h << ' ' << new_y_i[0] << '\n';
            y_i = new_y_i;
        }
        file.close();
    }

    void setFunction(std::vector<double (*)(double, std::vector<double>)> _functions) {
        functions = std::move(_functions);
    }

private:
    double x_max;
    double n;
    double h;
    double x_0;
    std::vector<double> y_0;
    std::vector<double (*)(double, std::vector<double>)> functions;
};

class RungeKutta5 {
public:
    RungeKutta5(double _x_max, double _n, double _x_0, std::vector<double> _y_0) : x_max(_x_max), n(_n), x_0(_x_0), y_0(std::move(_y_0)) {
        h = (x_max - _x_0) / n;
    }

    void printDataInFileFor123(std::vector<std::string> file_names) {
        std::vector<std::ofstream> files(file_names.size());
        for (int i = 0; i < file_names.size(); i++) {
            files[i].open(file_names[i]);
            if (!files[i]) {
                std::cout << "Error\n";
            }
        }
        std::vector<double> y_i = y_0;
        for (int i = 0; i < y_i.size(); i++) {
            files[i] << x_0 << ' ' << y_0[i] << '\n';
        }
        for (int i = 0; i < n; i++) {
            double x_i = x_0 + i * h;
            std::vector<double> new_y_i(y_i.size());
            std::vector<double> K(5);
            for (int j = 0; j < y_i.size(); j++) {
                K[0] = h * functions[j](x_i, y_i);

                std::vector<double> y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 3;
                K[1] = h * functions[j](x_i + h / 3, y_i_temp);

                y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 6 + K[1] / 6;
                K[2] = h * functions[j](x_i + h / 3, y_i_temp);

                y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 8 + (K[2] * 3) / 8;
                K[3] = h * functions[j](x_i + h / 2, y_i_temp);

                y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 2 - (3 * K[2]) / 2 + 2 * K[3];
                K[4] = h * functions[j](x_i + h, y_i_temp);

                new_y_i[j] = y_i[j] + (1 / 6.) * (K[0] + 4 * K[3] + K[4]);
                files[j] << x_i + h << ' ' << new_y_i[j] << '\n';
            }
            y_i = new_y_i;
        }
        for (int i = 0; i < y_i.size(); i++) {
            files[i].close();
        }
    }

    void printDataInFileFor4(const std::string& file_name) {
        std::ofstream file;
        file.open(file_name);
        if (!file) {
            std::cout << "Error\n";
        }
        std::vector<double> y_i = y_0;
        file << x_0 << ' ' << y_i[0] << '\n';
        for (int i = 0; i < n; i++) {
            double x_i = x_0 + i * h;
            std::vector<double> new_y_i(y_i.size());
            std::vector<double> K(5);
            for (int j = 0; j < y_i.size(); j++) {
                K[0] = h * functions[j](x_i, y_i);

                std::vector<double> y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 3;
                K[1] = h * functions[j](x_i + h / 3, y_i_temp);

                y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 6 + K[1] / 6;
                K[2] = h * functions[j](x_i + h / 3, y_i_temp);

                y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 8 + (K[2] * 3) / 8;
                K[3] = h * functions[j](x_i + h / 2, y_i_temp);

                y_i_temp = y_i;
                for (auto &k : y_i_temp) k += K[0] / 2 - (3 * K[2]) / 2 + 2 * K[3];
                K[4] = h * functions[j](x_i + h, y_i_temp);

                new_y_i[j] = y_i[j] + (1 / 6.) * (K[0] + 4 * K[3] + K[4]);
            }
            file << x_i + h << ' ' << new_y_i[0] << '\n';
            y_i = new_y_i;
        }
        file.close();
    }

    void setFunction(std::vector<double (*)(double, std::vector<double>)> _functions) {
        functions = std::move(_functions);
    }

private:
    double x_max;
    double n;
    double h;
    double x_0;
    std::vector<double> y_0;
    std::vector<double (*)(double, std::vector<double>)> functions;
};

class AdamsBashfort4 {
public:
    AdamsBashfort4(double _x_max, double _n, double _x_0, std::vector<double> _y_0) : x_max(_x_max), n(_n), x_0(_x_0), y_0(std::move(_y_0)) {
        h = (x_max - _x_0) / n;
    }

    void printDataInFileFor123(std::vector<std::string> file_names) {
        std::vector<std::ofstream> files(file_names.size());
        for (int i = 0; i < file_names.size(); i++) {
            files[i].open(file_names[i]);
            if (!files[i]) {
                std::cout << "Error\n";
            }
        }
        std::vector<std::vector<double>> Y(4, std::vector<double> (y_0.size()));
        std::vector<double> X(4);
        Y[0] = y_0;
        X[0] = x_0;
        for (int i = 0; i < 3; i++) {
            X[i + 1] = X[i] + h;
            std::vector<double> new_y_i(Y[i].size());
            std::vector<double> K(4);
            for (int j = 0; j < Y[i].size(); j++) {
                K[0] = h * functions[j](X[i], Y[i]);

                std::vector<double> y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[0] / 2;
                K[1] = h * functions[j](X[i] + h / 2, y_i_temp);

                y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[1] / 2;
                K[2] = h * functions[j](X[i] + h / 2, y_i_temp);

                y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[2];
                K[3] = h * functions[j](X[i] + h, y_i_temp);

                new_y_i[j] = Y[i][j] + (1 / 6.) * (K[0] + 2 * K[1] + 2 * K[2] + K[3]);
                files[j] << X[i] << ' ' << Y[i][j] << '\n';
            }
            Y[i + 1] = new_y_i;
        }
        for (int i = 0; i < Y[0].size(); i++) {
            files[i] << X[3] << ' ' << Y[3][i] << '\n';
        }
        for (int i = 4; i < n; i++) {
            std::vector<double> new_Y(Y[0].size());
            for (int j = 0; j < Y[0].size(); j++) {
                new_Y[j] = Y[3][j] + (h / 24) * (55 * functions[j](X[3], Y[3]) - 59 * functions[j](X[2], Y[2]) +
                        37 * functions[j](X[1], Y[1]) - 9 * functions[j](X[0], Y[0]));
                files[j] << X[3] + h << ' ' << new_Y[j] << '\n';
            }
            for (int k = 0; k < 3; k++) {
                X[k] = X[k + 1];
                Y[k] = Y[k + 1];
            }
            X[3] += h;
            Y[3] = new_Y;
        }
        for (auto & file : files) {
            file.close();
        }
    }

    void printDataInFileFor4(const std::string& file_name) {
        std::ofstream file;
        file.open(file_name);
        if (!file) {
            std::cout << "Error\n";
        }
        std::vector<std::vector<double>> Y(4, std::vector<double> (y_0.size()));
        std::vector<double> X(4);
        Y[0] = y_0;
        X[0] = x_0;
        for (int i = 0; i < 3; i++) {
            X[i + 1] = X[i] + h;
            std::vector<double> new_y_i(Y[i].size());
            std::vector<double> K(4);
            for (int j = 0; j < Y[i].size(); j++) {
                K[0] = h * functions[j](X[i], Y[i]);

                std::vector<double> y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[0] / 2;
                K[1] = h * functions[j](X[i] + h / 2, y_i_temp);

                y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[1] / 2;
                K[2] = h * functions[j](X[i] + h / 2, y_i_temp);

                y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[2];
                K[3] = h * functions[j](X[i] + h, y_i_temp);

                new_y_i[j] = Y[i][j] + (1 / 6.) * (K[0] + 2 * K[1] + 2 * K[2] + K[3]);
            }
            file << X[i] << ' ' << Y[i][0] << '\n';
            Y[i + 1] = new_y_i;
        }
        file << X[3] << ' ' << Y[3][0] << '\n';
        for (int i = 4; i < n; i++) {
            std::vector<double> new_Y(Y[0].size());
            for (int j = 0; j < Y[0].size(); j++) {
                new_Y[j] = Y[3][j] + (h / 24) * (55 * functions[j](X[3], Y[3]) - 59 * functions[j](X[2], Y[2]) +
                                                 37 * functions[j](X[1], Y[1]) - 9 * functions[j](X[0], Y[0]));
            }
            file << X[3] + h << ' ' << new_Y[0] << '\n';
            for (int k = 0; k < 3; k++) {
                X[k] = X[k + 1];
                Y[k] = Y[k + 1];
            }
            X[3] += h;
            Y[3] = new_Y;
        }
        file.close();
    }

    void setFunction(std::vector<double (*)(double, std::vector<double>)> _functions) {
        functions = std::move(_functions);
    }

private:
    double x_max;
    double n;
    double h;
    double x_0;
    std::vector<double> y_0;
    std::vector<double (*)(double, std::vector<double>)> functions;
};

class AdamsMoulton4 {
public:
    AdamsMoulton4(double _x_max, double _n, double _x_0, std::vector<double> _y_0) : x_max(_x_max), n(_n), x_0(_x_0), y_0(std::move(_y_0)) {
        h = (x_max - _x_0) / n;
        epsilon = 0.0000001;
    }

    void setEpsilon(double _epsilon) {
        epsilon = _epsilon;
    }

    void printDataInFileFor123(std::vector<std::string> file_names) {
        std::vector<std::ofstream> files(file_names.size());
        for (int i = 0; i < file_names.size(); i++) {
            files[i].open(file_names[i]);
            if (!files[i]) {
                std::cout << "Error\n";
            }
        }
        std::vector<std::vector<double>> Y(4, std::vector<double> (y_0.size()));
        std::vector<double> X(4);
        Y[0] = y_0;
        X[0] = x_0;
        for (int i = 0; i < 3; i++) {
            X[i + 1] = X[i] + h;
            std::vector<double> new_y_i(Y[i].size());
            std::vector<double> K(4);
            for (int j = 0; j < Y[i].size(); j++) {
                K[0] = h * functions[j](X[i], Y[i]);

                std::vector<double> y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[0] / 2;
                K[1] = h * functions[j](X[i] + h / 2, y_i_temp);

                y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[1] / 2;
                K[2] = h * functions[j](X[i] + h / 2, y_i_temp);

                y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[2];
                K[3] = h * functions[j](X[i] + h, y_i_temp);

                new_y_i[j] = Y[i][j] + (1 / 6.) * (K[0] + 2 * K[1] + 2 * K[2] + K[3]);
                files[j] << X[i] << ' ' << Y[i][j] << '\n';
            }
            Y[i + 1] = new_y_i;
        }
        for (int i = 0; i < Y[0].size(); i++) {
            files[i] << X[3] << ' ' << Y[3][i] << '\n';
        }
        for (int i = 4; i < n; i++) {
            std::vector<double> new_Y(Y[0].size());
            std::vector<double> next_Y(new_Y.size());
            for (int j = 0; j < Y[0].size(); j++) {
                new_Y[j] = Y[3][j] + (h / 24) * (55 * functions[j](X[3], Y[3]) - 59 * functions[j](X[2], Y[2]) +
                                                 37 * functions[j](X[1], Y[1]) - 9 * functions[j](X[0], Y[0]));
                do {
                    next_Y[j] =
                            Y[3][j] + (h / 24) * (9 * functions[j](X[3] + h, new_Y) + 19 * functions[j](X[3], Y[3]) -
                                                  5 * functions[j](X[2], Y[2]) + functions[j](X[1], Y[1]));
                    new_Y[j] = next_Y[j];
                } while (std::abs(next_Y[j] - new_Y[j]) > epsilon);
                files[j] << X[3] + h << ' ' << next_Y[j] << '\n';
            }
            for (int k = 0; k < 3; k++) {
                X[k] = X[k + 1];
                Y[k] = Y[k + 1];
            }
            X[3] += h;
            Y[3] = new_Y;
        }
        for (auto & file : files) {
            file.close();
        }
    }

    void printDataInFileFor4(const std::string& file_name) {
        std::ofstream file;
        file.open(file_name);
        if (!file) {
            std::cout << "Error\n";
        }
        std::vector<std::vector<double>> Y(4, std::vector<double> (y_0.size()));
        std::vector<double> X(4);
        Y[0] = y_0;
        X[0] = x_0;
        for (int i = 0; i < 3; i++) {
            X[i + 1] = X[i] + h;
            std::vector<double> new_y_i(Y[i].size());
            std::vector<double> K(4);
            for (int j = 0; j < Y[i].size(); j++) {
                K[0] = h * functions[j](X[i], Y[i]);

                std::vector<double> y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[0] / 2;
                K[1] = h * functions[j](X[i] + h / 2, y_i_temp);

                y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[1] / 2;
                K[2] = h * functions[j](X[i] + h / 2, y_i_temp);

                y_i_temp = Y[i];
                for (auto &k : y_i_temp) k += K[2];
                K[3] = h * functions[j](X[i] + h, y_i_temp);

                new_y_i[j] = Y[i][j] + (1 / 6.) * (K[0] + 2 * K[1] + 2 * K[2] + K[3]);
            }
            file << X[i] << ' ' << Y[i][0] << '\n';
            Y[i + 1] = new_y_i;
        }
        file << X[3] << ' ' << Y[3][0] << '\n';
        for (int i = 4; i < n; i++) {
            std::vector<double> new_Y(Y[0].size());
            std::vector<double> next_Y(new_Y.size());
            for (int j = 0; j < Y[0].size(); j++) {
                new_Y[j] = Y[3][j] + (h / 24) * (55 * functions[j](X[3], Y[3]) - 59 * functions[j](X[2], Y[2]) +
                                                 37 * functions[j](X[1], Y[1]) - 9 * functions[j](X[0], Y[0]));
                do {
                    next_Y[j] =
                            Y[3][j] + (h / 24) * (9 * functions[j](X[3] + h, new_Y) + 19 * functions[j](X[3], Y[3]) -
                                                  5 * functions[j](X[2], Y[2]) + functions[j](X[1], Y[1]));
                    new_Y[j] = next_Y[j];
                } while (std::abs(next_Y[j] - new_Y[j]) > epsilon);
            }
            file << X[3] + h << ' ' << new_Y[0] << '\n';
            for (int k = 0; k < 3; k++) {
                X[k] = X[k + 1];
                Y[k] = Y[k + 1];
            }
            X[3] += h;
            Y[3] = new_Y;
        }
        file.close();
    }

    void setFunction(std::vector<double (*)(double, std::vector<double>)> _functions) {
        functions = std::move(_functions);
    }

private:
    double epsilon;
    double x_max;
    double n;
    double h;
    double x_0;
    std::vector<double> y_0;
    std::vector<double (*)(double, std::vector<double>)> functions;
};

double function_first(double x, std::vector<double> y) {
    return y[0] * y[0] * sin(x * x + y[0]);
}

double function_second_1(double x, std::vector<double> y) {
    return pow(cos(x + y[1]), 2);
}

double function_second_2(double x, std::vector<double> y) {
    return sin(x - y[0]) / 2;
}

double function_third_1(double x, std::vector<double> y) {
    return log(pow(x, 2) + pow(y[0], 2));
}

double function_third_2(double x, std::vector<double> y) {
    return atan(x * y[0] * y[2]);
}

double function_third_3(double x, std::vector<double> y) {
    return sin(atan(y[0] * y[1]));
}

double function_fourth_1(double x, std::vector<double> y) {
    return y[1];
}

double function_fourth_2(double x, std::vector<double> y) {
    return y[2];
}
// Expressed function
double function_fourth_3(double x, std::vector<double> y) {
    //return cos(x / (x + 2)) - 2*y[2] - exp(x) * y[1] - atan(pow(x, 2) + 4) * y[0];
    //return exp(x)*cos(x)-pow(x,2)*y[1]-x*y[0];
    return pow(x, 2) * sin(x) - asinh(x) * y[2] - tan(pow(x, 2) + 4) * y[1] - sqrt(pow(x, 2) + 4) * y[0];
}

void euler() {
    std::string first_file_name = "euler_first_number";
    std::string second_file_1_name = "euler_second_number_1";
    std::string second_file_2_name = "euler_second_number_2";
    std::string third_file_1_name = "euler_third_number_1";
    std::string third_file_2_name = "euler_third_number_2";
    std::string third_file_3_name = "euler_third_number_3";
    std::string fourth_file_name = "euler_fourth_number";

    // Euler. Task 1.
    GnuplotPipe gp1;
    Euler euler_first(4, 100000, 0, {0.5});
    euler_first.setFunction({function_first});
    euler_first.printDataInFileFor123({first_file_name});
    gp1.sendLine("plot \"" + first_plot_file_path + first_file_name + R"(" with line title "y1")");
    //*****************************************************

    // Euler. Task 2.
    GnuplotPipe gp2;
    Euler euler_second(3, 100000, 0, {0.3, -0.6});
    euler_second.setFunction({function_second_1, function_second_2});
    euler_second.printDataInFileFor123({second_file_1_name, second_file_2_name});
    gp2.sendLine("plot \"" + first_plot_file_path + second_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + second_file_2_name + R"(" with lines title "y2")");
    //*****************************************************

    // Euler. Task 3.
    GnuplotPipe gp3;
    Euler euler_third(3, 100000, -3, {-1, -7, 2});
    euler_third.setFunction({function_third_1, function_third_2, function_third_3});
    euler_third.printDataInFileFor123({third_file_1_name, third_file_2_name, third_file_3_name});
    gp3.sendLine("plot \"" + first_plot_file_path + third_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + third_file_2_name + R"(" with lines title "y2", ")"
                 + first_plot_file_path + third_file_3_name + R"(" with lines title "y3")");
    //*****************************************************

    // Euler. Task 4.
    GnuplotPipe gp4;
    Euler euler_fourth(0.5, 100000, -0.5, {2, -2, -5});
    euler_fourth.setFunction({function_fourth_1, function_fourth_2, function_fourth_3});
    euler_fourth.printDataInFileFor4({fourth_file_name});
    gp4.sendLine("plot \"" + first_plot_file_path + fourth_file_name + R"(" with lines title "y1")");
    //*****************************************************
}

void runge_kutta_4() {
    std::string first_file_name = "runge_kutta_4_first_number";
    std::string second_file_1_name = "runge_kutta_4_second_number_1";
    std::string second_file_2_name = "runge_kutta_4_second_number_2";
    std::string third_file_1_name = "runge_kutta_4_third_number_1";
    std::string third_file_2_name = "runge_kutta_4_third_number_2";
    std::string third_file_3_name = "runge_kutta_4_third_number_3";
    std::string fourth_file_name = "runge_kutta_4_fourth_number";

    // Runge Kutta 4. Task 1.
    GnuplotPipe gp1;
    RungeKutta4 euler_first(4, 100000, 0, {0.5});
    euler_first.setFunction({function_first});
    euler_first.printDataInFileFor123({first_file_name});
    gp1.sendLine("plot \"" + first_plot_file_path + first_file_name + R"(" with line title "y1")");
    //*****************************************************

    // Runge Kutta 4. Task 2.
    GnuplotPipe gp2;
    RungeKutta4 euler_second(3, 100000, 0, {0.3, -0.6});
    euler_second.setFunction({function_second_1, function_second_2});
    euler_second.printDataInFileFor123({second_file_1_name, second_file_2_name});
    gp2.sendLine("plot \"" + first_plot_file_path + second_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + second_file_2_name + R"(" with lines title "y2")");
    //*****************************************************

    // Runge Kutta 4. Task 3.
    GnuplotPipe gp3;
    RungeKutta4 euler_third(3, 100000, -3, {-1, -7, 2});
    euler_third.setFunction({function_third_1, function_third_2, function_third_3});
    euler_third.printDataInFileFor123({third_file_1_name, third_file_2_name, third_file_3_name});
    gp3.sendLine("plot \"" + first_plot_file_path + third_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + third_file_2_name + R"(" with lines title "y2", ")"
                 + first_plot_file_path + third_file_3_name + R"(" with lines title "y3")");
    //*****************************************************

    // Runge Kutta 4. Task 4.
    GnuplotPipe gp4;
    RungeKutta4 euler_fourth(0.5, 100000, -0.5, {2, -2, -5});
    euler_fourth.setFunction({function_fourth_1, function_fourth_2, function_fourth_3});
    euler_fourth.printDataInFileFor4({fourth_file_name});
    gp4.sendLine("plot \"" + first_plot_file_path + fourth_file_name + R"(" with lines title "y1")");
    //*****************************************************
}

void runge_kutta_5() {
    std::string first_file_name = "runge_kutta_5_first_number";
    std::string second_file_1_name = "runge_kutta_5_second_number_1";
    std::string second_file_2_name = "runge_kutta_5_second_number_2";
    std::string third_file_1_name = "runge_kutta_5_third_number_1";
    std::string third_file_2_name = "runge_kutta_5_third_number_2";
    std::string third_file_3_name = "runge_kutta_5_third_number_3";
    std::string fourth_file_name = "runge_kutta_5_fourth_number";

    // Runge Kutta 5. Task 1.
    GnuplotPipe gp1;
    RungeKutta5 euler_first(4, 100000, 0, {0.5});
    euler_first.setFunction({function_first});
    euler_first.printDataInFileFor123({first_file_name});
    gp1.sendLine("plot \"" + first_plot_file_path + first_file_name + R"(" with line title "y1")");
    //*****************************************************

    // Runge Kutta 5. Task 2.
    GnuplotPipe gp2;
    RungeKutta5 euler_second(3, 100000, 0, {0.3, -0.6});
    euler_second.setFunction({function_second_1, function_second_2});
    euler_second.printDataInFileFor123({second_file_1_name, second_file_2_name});
    gp2.sendLine("plot \"" + first_plot_file_path + second_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + second_file_2_name + R"(" with lines title "y2")");
    //*****************************************************

    // Runge Kutta 5. Task 3.
    GnuplotPipe gp3;
    RungeKutta5 euler_third(3, 100000, -3, {-1, -7, 2});
    euler_third.setFunction({function_third_1, function_third_2, function_third_3});
    euler_third.printDataInFileFor123({third_file_1_name, third_file_2_name, third_file_3_name});
    gp3.sendLine("plot \"" + first_plot_file_path + third_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + third_file_2_name + R"(" with lines title "y2", ")"
                 + first_plot_file_path + third_file_3_name + R"(" with lines title "y3")");
    //*****************************************************

    // Runge Kutta 5. Task 4.
    GnuplotPipe gp4;
    RungeKutta5 euler_fourth(0.5, 100000, -0.5, {2, -2, -5});
    euler_fourth.setFunction({function_fourth_1, function_fourth_2, function_fourth_3});
    euler_fourth.printDataInFileFor4({fourth_file_name});
    gp4.sendLine("plot \"" + first_plot_file_path + fourth_file_name + R"(" with lines title "y1")");
    //*****************************************************
}

void adams_bashfort_4() {
    std::string first_file_name = "adams_bashfort_4_first_number";
    std::string second_file_1_name = "adams_bashfort_4_second_number_1";
    std::string second_file_2_name = "adams_bashfort_4_second_number_2";
    std::string third_file_1_name = "adams_bashfort_4_third_number_1";
    std::string third_file_2_name = "adams_bashfort_4_third_number_2";
    std::string third_file_3_name = "adams_bashfort_4_third_number_3";
    std::string fourth_file_name = "adams_bashfort_4_fourth_number";

    // Adams Bashfort 4. Task 1.
    GnuplotPipe gp1;
    AdamsBashfort4 euler_first(4, 100000, 0, {0.5});
    euler_first.setFunction({function_first});
    euler_first.printDataInFileFor123({first_file_name});
    gp1.sendLine("plot \"" + first_plot_file_path + first_file_name + R"(" with line title "y1")");
    //*****************************************************

    // Adams Bashfort 4. Task 2.
    GnuplotPipe gp2;
    AdamsBashfort4 euler_second(3, 100000, 0, {0.3, -0.6});
    euler_second.setFunction({function_second_1, function_second_2});
    euler_second.printDataInFileFor123({second_file_1_name, second_file_2_name});
    gp2.sendLine("plot \"" + first_plot_file_path + second_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + second_file_2_name + R"(" with lines title "y2")");
    //*****************************************************

    // Adams Bashfort 4. Task 3.
    GnuplotPipe gp3;
    AdamsBashfort4 euler_third(3, 100000, -3, {-1, -7, 2});
    euler_third.setFunction({function_third_1, function_third_2, function_third_3});
    euler_third.printDataInFileFor123({third_file_1_name, third_file_2_name, third_file_3_name});
    gp3.sendLine("plot \"" + first_plot_file_path + third_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + third_file_2_name + R"(" with lines title "y2", ")"
                 + first_plot_file_path + third_file_3_name + R"(" with lines title "y3")");
    //*****************************************************

    // Adams Bashfort 4. Task 4.
    GnuplotPipe gp4;
    AdamsBashfort4 euler_fourth(0.5, 100000, -0.5, {2, -2, -5});
    euler_fourth.setFunction({function_fourth_1, function_fourth_2, function_fourth_3});
    euler_fourth.printDataInFileFor4({fourth_file_name});
    gp4.sendLine("plot \"" + first_plot_file_path + fourth_file_name + R"(" with lines title "y1")");
    //*****************************************************
}

void adams_moulton_4() {
    std::string first_file_name = "adams_moulton_4_first_number";
    std::string second_file_1_name = "adams_moulton_4_second_number_1";
    std::string second_file_2_name = "adams_moulton_4_second_number_2";
    std::string third_file_1_name = "adams_moulton_4_third_number_1";
    std::string third_file_2_name = "adams_moulton_4_third_number_2";
    std::string third_file_3_name = "adams_moulton_4_third_number_3";
    std::string fourth_file_name = "adams_moulton_4_fourth_number";

    // Adams Moulton 4. Task 1.
    GnuplotPipe gp1;
    AdamsMoulton4 euler_first(4, 100000, 0, {0.5});
    euler_first.setFunction({function_first});
    euler_first.printDataInFileFor123({first_file_name});
    gp1.sendLine("plot \"" + first_plot_file_path + first_file_name + R"(" with line title "y1")");
    //*****************************************************

    // Adams Moulton 4. Task 2.
    GnuplotPipe gp2;
    AdamsMoulton4 euler_second(3, 100000, 0, {0.3, -0.6});
    euler_second.setFunction({function_second_1, function_second_2});
    euler_second.printDataInFileFor123({second_file_1_name, second_file_2_name});
    gp2.sendLine("plot \"" + first_plot_file_path + second_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + second_file_2_name + R"(" with lines title "y2")");
    //*****************************************************

    // Adams Moulton 4. Task 3.
    GnuplotPipe gp3;
    AdamsMoulton4 euler_third(3, 100000, -3, {-1, -7, 2});
    euler_third.setFunction({function_third_1, function_third_2, function_third_3});
    euler_third.printDataInFileFor123({third_file_1_name, third_file_2_name, third_file_3_name});
    gp3.sendLine("plot \"" + first_plot_file_path + third_file_1_name + R"(" with lines title "y1", ")"
                 + first_plot_file_path + third_file_2_name + R"(" with lines title "y2", ")"
                 + first_plot_file_path + third_file_3_name + R"(" with lines title "y3")");
    //*****************************************************

    // Adams Moulton 4. Task 4.
    GnuplotPipe gp4;
    AdamsMoulton4 euler_fourth(0.5, 100000, -0.5, {2, -2, -5});
    euler_fourth.setFunction({function_fourth_1, function_fourth_2, function_fourth_3});
    euler_fourth.printDataInFileFor4({fourth_file_name});
    gp4.sendLine("plot \"" + first_plot_file_path + fourth_file_name + R"(" with lines title "y1")");
    //*****************************************************
}

int main(){
    euler();

    return 0;
}