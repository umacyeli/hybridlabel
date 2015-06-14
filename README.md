# hybridlabel
Hybrid labeling(HYBRID) is our original labeling technique for shortest-path and distance queries on a variety of graphs. HYBRID captures characteristics from several existing state-of-the-art techniques such as CH, PLL and HopDoubling and it provides with a controllable performance tuning which can efficiently trade off between fast query and  low space consumption.

# description
decompose() identifies two categories labels(structual+shortcuts vs. guidences), lightenLabel() consistently removes a single label from the label lists of each vertices.

In lightenLabel(), there is another method toTree() which can tranfer the original PLL labels into a shortest path trees. popTree() pops a removed label for each tree.

Each time we remove 1 label from each lists and check the decreased efficiency of query time.


#We need to check(let's call this a TODO list):

##Dropping Stage
Uniform distribution over vertices(every vertices lose the same number of label unless one has nothing to lose)
1) "Top Down";
2) "Bottom Up";

skew distribution(what is the goal function? pop the least shortest-path coverage label, 5%-10%-....-80% vs. query performance)
3) heustic(greedy)

##Query Stage
Any ideas to boost up the performance?

Currently(query u,v):
prequites:
1. id-sorted remaining labels(nlabel, nlabeld, glabel, glabeld);
2. distance-sorted expanding parents(eparent, eparentd); 

procedures:
1. Merge join the remaining lables, to have a MIN(not really MIN);
2. Check the smallest expanding parents, if MIN is not larger than both smallest parents, return MIN, or else push all expanding parents into respective heaps;
3. Bidirectional dijkstra...

Maybe optimization:
1. A\*?;
2. Landmark?

##Comprehensive experiments
1. greedy order 2-hop(Cohen);
2. TNR?;
3. Space overhead scale(our vs. ch. vs pll) expecting: pll going widely(quadratic) while our is consistent with ch(linear)
4. Query time scale(ours vs. ch vs. pll) expecting; ch going widely(quadratic) while our is consistent with pll(linear)


##Some justification
###Why we insist on degree ordering in processing?
For example, brightkite
1. the performance gap between cohen's method and pll;
2. the performance gap between CH and our minimal labels.(0.8m vs. 1m shortcuts)
3. Also, how about AH(although it only works on road networks...)

Are these two gaps highly correlated? I think so...
