#ifndef CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_
#define CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>
#include <stdexcept>

constexpr double EPS = 1e-7;

class S21Matrix {
 public:
  // Constructors & Destructor
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();
  // Get & Set methods
  int GetRows() const;
  void SetRows(int rows);
  int GetCols() const;
  void SetCols(int cols);

  // Member functions
  void SumMatrix(const S21Matrix& other);
  bool EqMatrix(const S21Matrix& other) const;
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  double Determinant() const;
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  S21Matrix InverseMatrix() const;

  // Operators
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(double x) const;
  bool operator==(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(double x);
  double& operator()(int rows, int cols);
  double operator()(int rows, int cols) const;

 private:
  // Attributes
  int rows_, cols_;
  double** matrix_;
  // Sub functions
  void CreateMatrix();
  void DeleteMatrix();
  void CopyMatrix(const S21Matrix& other);
  void MatrixMinors(const S21Matrix& other, int rows_i, int cols_j) const;
  double MinorsDeterminants() const;
};

#endif  // CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_
