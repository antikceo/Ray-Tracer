
Note: Some initial board states take a long time to search, so if it looks like nothing is happening, it's just searching. If it takes too long, just reload, maybe the random case is a really hard case.

Program outputs results for randomly generated states as follows:



To test the program, enter your input in the form:
	(8-puzzle '(3 1 2 6 4 5 0 7 8) 'manhattan) or
	(8-puzzle '(3 1 2 6 4 5 0 7 8) 'misplaced) or
	(8-puzzle '(3 1 2 6 4 5 0 7 8) 'extracredit)

The goal state is '(0 1 2 3 4 5 6 7 8)

Manhattan functions, compute the manhattan distance by summing the row and column distances of each tile from the goal.

Misplaced checks each tile in the state and if the position does not match the goal position, 1 is added to the heuristic.

Concept for extra credit gotten from web.mit.edu/6.034/wwwbob/EightPuzzle.pdf

The idea is to create a tile reversal penalty will add one penalty move for a direct  reversal of tiles i.e. the board:

213
804
765
would have a reversal penalty of 2, because 2 is reversed with 1 and 1 is reversed with two. Reversal penalties will always come in pairs.

Reversed tiles are much more difficult to deal with (they must go around each other), for every two reversed tiles, this will require a minimum of 2 moves.

On some runs, the Manhattan actually did slightly better than the penalty heuristic, but the penalty heuristic dominated in most cases.

To run random case: (random-case)

Sample run output: 

Initial board 1..
5 3 7
4 0 8
6 2 1
 
  Manhattan: 2553
  Misplaced: 46808
  ExtraCredit: 2312

Initial board 2..
3 5 2
1 4 6
7 8 0
 
  Manhattan: 1593
  Misplaced: 9577
  ExtraCredit: 1120

Initial board 3..
5 7 3
1 0 8
6 2 4
 
  Manhattan: 806
  Misplaced: 11768
  ExtraCredit: 573

Initial board 4..
1 7 6
8 3 4
0 5 2
 
  Manhattan: 1006
  Misplaced: 23830
  ExtraCredit: 861

Initial board 5..
0 6 1
8 3 7
4 5 2
 
  Manhattan: 679
  Misplaced: 8403
  ExtraCredit: 301
