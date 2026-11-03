import sys

sys.path.append('../inputs')
sys.path.append('../solvers')
sys.path.append('../kinetics')
sys.path.append('../misc')
sys.path.append('.')

class linear_ldf_c :

    def Linear(self, q,s):
        K2 = s.Expand(s.K,s.Nz,1)
        %Nz x Nc
        c = q./(s.qm*K2)
        return c

class linear_ldf_q :

    def Linear(self, c,s):
        K2 = s.Expand(s.K,s.Nz,1); %Nz x Nc
        q = c*(s.qm*K2);
        return q

class linear_pore_c :

    def Linear(c,q,s):
        c = q/0.000001
        return c

class linear_pore_q : 

    def Linear(c,s):
        q = 0*c;
        return q

