## Lab 2
### Creating a DoFHandler and visualising sparsity patterns

author: Luca Heltai (luca.heltai@sissa.it)

Some useful resources
https://www.dealii.org/current/doxygen/deal.II/step_2.html https://www.dealii.org/current/doxygen/deal.II/namespaceDoFRenumbering.html

### Using step-2 as a base:
1. Compile and run this tutorial, and inspect at the output.
2. Look at the generated sparsity patterns (Firefox, for example).
2. Investigate:
	- How does the pattern change if you increase the polynomial degree from 1 to 2, or to 3?
	- How does the pattern change if you use a globally refined (say 3 times) unit square?
	- Are these patterns symmetric? Why/why not?
	- How many entries per row in the sparsity pattern do you expect for a Q1 element (assuming four cells are around each vertex)?
	- Check that this is true for the mesh in (b) (look for row_length(i) and output them for each row).
	- Can you construct a 2d mesh (without hanging nodes) that has a row with more entries?
	- How many entries per row in the sparsity pattern are there for Q2 and Q3 elements, again assuming four cells around each vertex?
	- Print all entries for row 42 for the original renumbered sparsity pattern.
	- Renumber the DoFs using the `boost::king_ordering` algorithm. What does the sparsity pattern look like now?
3. Additional tasks
	- Compute and output statistics like the number of unknowns, bandwidth of the sparsity pattern, average number of entries per row, and fill ratio.
	- Investigate the other appropriate DoF renumbering schemes. Which one produces the most banded structure?
	- Repeat the above for increasing refinement levels. Which is the most effcient scheme (lets say, in terms of bandwidth reduction versus computational time expended)? You can get an estimate of the time for this operation like this:
`$ time ./step-2`
	- What happens if you change the mesh from 2d to 3d?
	- Investigate the sparsity patterns generated for other types of `FiniteElements` of varying order.
