--- Exam Project 14 ---
by William Jedig

My evaluation for my score on this exam project is (10/10), as I feel I have completed all the tasks to a high standard.

This project involves computing the Cholesky Decomposition of real symmetric positive-definite matrix, which I will call matrix A.

I implemented the Cholesky–Banachiewicz algorithm, in accordance with the Wikipedia article.

I then implemented a simple method for generating real symmetric positive definite matrices, by generating a random square matrix M with values between [0,1]. I then generated matrix A by computing the product M*(M^T), which produces a real symmetric positive definite matrix per: 

https://www.physicsforums.com/threads/all-the-ways-to-build-positive-definite-matrices.561438/

To implement the linear equation solver, I implemented the following algorithm from the Wikipedia article:

Goal: Solve Ax=b
1. Decompose A to A = LL^T
2. Solve for y in Ly=b, via forward substitution
3. Solve for x in (L^T)x=y via back substitution
4. return x

To implement the determinant finder, I used the fact that the determinant of a triangular matrix is the product of the diagonals. Since the diagonals of L and L^T are identical, this problem reduces to:

Det(A) = Det(LL^T) = Det(L)Det(L^T) = (Det(L))^2,

i.e. multiplication of the diagonal elements, followed by squaring the result.

Finally, to implement inverse solver, I used the fact that the inverse is given by the matrix of columns xi in:

A*xi = ei

Thus this problem reduces to simply using the linear equation solver, on the set of unit vector ei corresponding to the dimension of the problem.
