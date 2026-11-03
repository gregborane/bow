import sys

sys.path.append(".")
sys.path.append("../inputs/")
sys.path.append("../kinetics")
sys.path.append("../misc")
sys.path.append("../solvers")

class langmuirLDF_c :

    def Langmuir(self, c,q,s)

        K2 = s.Expand(s.K,s.Nz,1); %Nz x Nc
        % qm2 = s.Expand(s.qm,s.Nz,1); %Nz x Nc
        sumi = s.Expand(sum(K2.*c,2),1,s.Nc);
        sumj = sumi - K2.*c;
        c = q.*(1+sumj)./(K2.*(s.qm-q));
        return c

class langmuirLDF_q :

    def Langmuir(self, c,s)
        K2 = s.Expand(s.K,s.Nz,1); %Nz x Nc
        kinf2D = s.Expand(s.kinf,s.Nz,1);
        sumi = s.Expand(sum(K2.*c,2),1,s.Nc);
        q = s.qm*K2.*c./(1+sumi) + kinf2D/s.theta.*c;
        return q

class langmuirPore_c :

    def Langmuir(self, c,q,s)
        #c is needed since ceq is dependence on other component c concentrations
        K3 = s.Expand(permute(s.K,[1 3 2]),s.Nz,s.Nr,1); %Nz x Nr x Nc
        qm3 = s.Expand(permute(s.qm,[1 3 2]),s.Nz,s.Nr,1); %Nz x Nr x Nc
        sumi = repmat(sum(K3.*c,3),1,1,s.Nc);
        sumj = sumi - K3.*c;
        c = q.*(1+sumj)./(K3.*(qm3-q));
        return c


class langmuirPore_q :
    
    def Langmuir(self, c,s):
        K3 = s.Expand(s.K,s.Nz,s.Nr,1); %Nz x Nr x Nc
        sumi = s.Expand(sum(K3.*c,3),1,1,s.Nc);
        q = s.qm.*K3.*c./(1+sumi);        
        return q
