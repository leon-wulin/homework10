HOMEWORK 9: DISTANCE FIELDS & FAST MARCHING METHOD


NAME:  <Lin Wu >


COLLABORATORS AND OTHER RESOURCES:
List the names of everyone you talked to about this assignment
(classmates, TAs, ALAC tutors, upperclassmen, students/instructor via
LMS, etc.), and all of the resources (books, online reference
material, etc.) you consulted in completing this assignment.

< w3cschool >

Remember: Your implementation for this assignment must be done on your
own, as described in "Academic Integrity for Homework" handout.



ESTIMATE OF # OF HOURS SPENT ON THIS ASSIGNMENT:  < 8 hours >



NAIVE ALGORITHM:

Order Notation:O(N^4) for a NxN image

Timing Experiment Data: 77289ms for 300*300

Discussion:the naive algorithm would have poor performance for large images.




IMPROVED ALGORITHM:

Order Notation: O(N^2 * K) for a NxN image, where K is the number of iterations until convergence.

Timing Experiment Data: 869ms for 300*300    58530ms for 1000*1000

Discussion:The improved algorithm iteratively updates the distance values until no changes are made. 



FAST MARCHING METHOD (with a map):

Order Notation: O(N^2 * log(N^2)) for a NxN image

Timing Experiment Data: 1408ms for 300*300  16054ms for 1000*1000

Discussion:The Fast Marching Method uses a priority queue to process pixels in order of increasing distance. 
This ensures that pixels are updated in a more optimal manner, reducing the number of iterations needed to converge. 	


DISTANCE FIELD VISUALIZATIONS FOR EXTRA CREDIT:
I am trying to use HSV format instead of RGB format for visualization because HSV format is more 
in line with human perception of color and color transitions are more natural. 
In theory, the appearance will also be better. And although I also used rainbow colors for visualization, 
I flipped the colors from warm to cool. However, based on the results, it may not be very successful.
I think it's not as good as the provided visualization function. I must admit that my artistic sense is not very good.
You can enter 'gradient' in the visualization method location to call the function I wrote.



FAST MARCHING METHOD (with a hash table) FOR EXTRA CREDIT:

Order Notation:

Timing Experiment Data:

Discussion:



MISC. COMMENTS TO GRADER:  
Optional, please be concise!






