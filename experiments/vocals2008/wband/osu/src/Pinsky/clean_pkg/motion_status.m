% Kongsberg motion status read from Ken and Sergio's spreadsheet and keyed into matlab.
% SPdeS 12/20/2011

statusday=(315:337)';
statushour=0:23;
Status=[[1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	3	3	3	3	3	3	4	4	0
3	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	3	3	0	0
3	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	4	4	4	0	4
0	0	0	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	3	3	3	3	3
3	3	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	4	0	0	0	0	0	0	4
4	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2
2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	0	0	0	0	0	0	4
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	4	4	4	4	4	4	4	4	0	0	0
0	4	0	4	4	4	4	4	4	0	4	0	2	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	4	4	4	4	4	4	4	4	4	4	0	0	4	3	4
3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	4	4	0	0	0	0	0	0	0	0	0	0
0	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1]];
% 0 OK
% 1 no data
% 2 bias (<1 degree)
% 3 noisy
% 4 motion adj. off or failed