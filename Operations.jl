E = sum(HH)
dv = E/NN
dc = E/MM
S = 4
K = 2S + 7
# RBP

C2V_up = (K+5)dc - 6

residuals = 4

residual_decay = 2

decaying = 3

real_value_comps = 2MM + 2dc + dc*dv - 4

V2C_up = 2dv

C2V_up_CR = 2


RBP_ops = C2V_up + 
        (dv-1)*V2C_up + 
        (dv-1)*(dc-1)*(C2V_up + residuals) +
        real_value_comps

RD_RBP_ops = C2V_up + decaying + 
             (dv-1)*V2C_up + 
             (dv-1)*(dc-1)*(C2V_up + residuals + residual_decay) + 
             real_value_comps

S1 = 16
S2 = 2

List_comps = 2*S2*((dv-1)*(dc-1) + S1)
update_list_ops = 6(S1-1) + 2(dv-1)*(dc-1)*(S1-1)

List_RBP_ops = C2V_up + decaying + 
             (dv-1)*V2C_up + 
             (dv-1)*(dc-1)*(C2V_up + residuals + residual_decay) + 
             List_comps +
             update_list_ops

CR_RBP_ops = dv*(C2V_up_CR + decaying) + 
             dv*V2C_up + 
             dv*(dc-1)*(C2V_up + residuals + residual_decay + 1) + 
             real_value_comps + dc

CDR_RBP_ops = dv*(C2V_up_CR + decaying) + 
              (dv-1)*V2C_up +
              (dv-1)*(dc-1)*(C2V_up + residuals + residual_decay) + 
              real_value_comps

E*RBP_ops
E*RD_RBP_ops
E*List_RBP_ops
(E-NN)*CR_RBP_ops
E*CDR_RBP_ops



display(1000*9*E*RD_RBP_ops/(1.2e9))
display(1000*11*E*RBP_ops/(1.2e9))
display(1000*7*E*RBP_ops/(1.2e9))
display(1000*6*E*List_RBP_ops/(1.2e9))
display(1000*4*(E-NN)*CR_RBP_ops/(1.2e9))
MEDIA = (E*CDR_RBP_ops + (E-NN)*CR_RBP_ops)/2
display(1000*4*MEDIA/(1.2e9))

display(1000*3*E*RBP_ops/(1.2e9))
display(1000*3*E*RD_RBP_ops/(1.2e9))
display(1000*5*E*RBP_ops/(1.2e9))
display(1000*2*E*RBP_ops/(1.2e9))
display(1000*2*E*List_RBP_ops/(1.2e9))
display(1000*2*(E-NN)*CR_RBP_ops/(1.2e9))
display(1000*2*E*CDR_RBP_ops*MEDIA/(1.2e9))

((3E + 2NN + MM)*8 + MM*NN + MM + NN)/1024

((2E + 2NN + S1 + S2)*8 + MM*NN + MM + NN)/1024

((4E + 2NN + MM)*8 + MM*NN + MM + NN)/1024

