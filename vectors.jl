v1 = [1;9;map(y->19+y,0:5)]#;map(y->12+y,0:4)]#;map(y->23+y,0:8)]
v2 = [0.5;2.0;map(y->2.0-(1.53/6)*y,1:6)]#;map(y->1.0+(0.4/5)*y,1:5)]
run_param_scen_cal(0.075,0,[3],[1],999,999,v2,v1,71,1.5,(1.5*1.3),1,[999],999,1,250)

v1 = [1;map(y->106+y,0:12)] ##72 112 0,3
v2 = [0.695;map(k->0.695-k*(0.28/13),1:13)]
run_param_scen_cal(0.075,0.67,[1],[999],999,999,v2,v1,71,1.5,(1.5*1.3),2,[82],999,1,250)

v1 = [1;map(y->169+y,0:22)]#;map(y->168+y,0:19)]
v2 = [1.0;map(y->1.0-(0.6/23)*y,1:23)]#;map(y->1.0-(0.599/20)*y,1:20)]
run_param_scen_cal(0.075,1.0,[1],[999],999,999,v2,v1,71,1.5,(1.5*1.3),3,[999],161,2,250)


#######################################
################ without distancing ###
########################################
v1 = [1;9;23]#;map(y->12+y,0:4)]#;map(y->23+y,0:8)]
v2 = [0.5;2.0;1.0]#;map(y->1.0+(0.4/5)*y,1:5)]
run_param_scen_cal(0.075,0,[3],[1],23,999,v2,v1,71,1.5,(1.5*1.3),1,[999],999,1,250)

v1 = [1;107]#;map(y->107+y,0:11)] ##72 112 0,3
v2 = [0.695;1.0]#;map(k->0.695-k*(0.28/12),1:12)]
run_param_scen_cal(0.075,0.67,[1],[999],107,999,v2,v1,71,1.5,(1.5*1.3),2,[82],999,1,250)

v1 = [1]#;map(y->168+y,0:19)]#;map(y->164+y,0:30)]
v2 = [1.0]#;map(y->1.0-(0.599/20)*y,1:20)]#;map(y->0.8-(0.64/31)*y,1:31)]
run_param_scen_cal(0.075,1.0,[1],[999],168,999,v2,v1,71,1.5,(1.5*1.3),3,[999],161,2,250)

#####################################################################3
############# No vaccination ###################################
#######################################################################
v1 = [1;9;map(y->19+y,0:5)]#;map(y->12+y,0:4)]#;map(y->23+y,0:8)]
v2 = [0.5;2.0;map(y->2.0-(1.53/6)*y,1:6)]#;map(y->1.0+(0.4/5)*y,1:5)]
run_param_scen_cal(0.075,0,[3],[1],999,999,v2,v1,999,1.5,(1.5*1.3),1,[999],999,1,250)

v1 = [1;map(y->107+y,0:11)] ##72 112 0,3
v2 = [0.695;map(k->0.695-k*(0.28/12),1:12)]
run_param_scen_cal(0.075,0.67,[1],[999],999,999,v2,v1,999,1.5,(1.5*1.3),2,[82],999,1,250)

v1 = [1;map(y->167+y,0:20)]#;map(y->168+y,0:19)]
v2 = [1.0;map(y->1.0-(0.60/21)*y,1:21)]#;map(y->1.0-(0.599/20)*y,1:20)]
run_param_scen_cal(0.075,1.0,[1],[999],999,999,v2,v1,999,1.5,(1.5*1.3),3,[999],161,2,250)


#####################################################################################
#####################################3 FAST VACCINATION ################################
#######################################################################################3
v1 = [1;9;map(y->19+y,0:5)]#;map(y->12+y,0:4)]#;map(y->23+y,0:8)]
v2 = [0.5;2.0;map(y->2.0-(1.53/6)*y,1:6)]#;map(y->1.0+(0.4/5)*y,1:5)]
run_param_scen_cal(0.075,0,[3],[1],999,999,v2,v1,71,1.5,(1.5*1.3),1,[999],999,1,250,"fast")

v1 = [1;map(y->107+y,0:11)] ##72 112 0,3
v2 = [0.695;map(k->0.695-k*(0.28/12),1:12)]
run_param_scen_cal(0.075,0.67,[1],[999],999,999,v2,v1,71,1.5,(1.5*1.3),2,[82],999,1,250,"fast")

v1 = [1;map(y->167+y,0:20)]#;map(y->168+y,0:19)]
v2 = [1.0;map(y->1.0-(0.60/21)*y,1:21)]#;map(y->1.0-(0.599/20)*y,1:20)]
run_param_scen_cal(0.075,1.0,[1],[999],999,999,v2,v1,71,1.5,(1.5*1.3),3,[999],161,2,250,"fast")



#######################################
################ without distancing and fast vaccination ###
########################################
v1 = [1;9;23]#;map(y->12+y,0:4)]#;map(y->23+y,0:8)]
v2 = [0.5;2.0;1.0]#;map(y->1.0+(0.4/5)*y,1:5)]
run_param_scen_cal(0.075,0,[3],[1],23,999,v2,v1,71,1.5,(1.5*1.3),1,[999],999,1,250,"fast")

v1 = [1;107]#;map(y->107+y,0:11)] ##72 112 0,3
v2 = [0.695;1.0]#;map(k->0.695-k*(0.28/12),1:12)]
run_param_scen_cal(0.075,0.67,[1],[999],107,999,v2,v1,71,1.5,(1.5*1.3),2,[82],999,1,250,"fast")

v1 = [1]#;map(y->168+y,0:19)]#;map(y->164+y,0:30)]
v2 = [1.0]#;map(y->1.0-(0.599/20)*y,1:20)]#;map(y->0.8-(0.64/31)*y,1:31)]
run_param_scen_cal(0.075,1.0,[1],[999],168,999,v2,v1,71,1.5,(1.5*1.3),3,[999],161,2,250,"fast")

