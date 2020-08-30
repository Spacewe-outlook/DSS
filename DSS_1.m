function [Gbest,Fit_min,fit_save,dim_save]= PHA_1_func(fhd,Dimension,Dim_Continuous,Max_Gen,VRmin,VRmax,varargin)
Nd = Dimension;
Max_gen = Max_Gen;

if Nd==1
    LB=VRmin;
    UB=VRmax;
else
    LB = VRmin*ones(1,Nd);
    UB = VRmax*ones(1,Nd);
end

Step_Coefficient = 100; 
Step_len_exp = Step_Coefficient;
% ������ά
Dim_Continuous = Dim_Continuous; %ÿ�ε�����������ֵ
% ��¼��ǰ�ƶ����ڼ�ά
Now_Dim_Slide_A = 1;
% ��¼״̬2�ƶ����ڼ�ά
Now_Slide_2 = 1;
% ����̶���ʼλ��
center = (LB+UB)/2;
A = (LB+center)/2;
% ������Ӧ��ֵ
fit_save=zeros(1,Max_gen);
%����ά�ȱ仯
dim_save =zeros(1,Max_gen);
% ����ֲ�����
local_history = [];
local_history = [local_history;A;];
% ȫ�����Ž�
Fit_min = Inf;
% ״̬��־
SF = 1;
temp_local_center = A;
s_num = 1;
%ά�ȱ仯��������
Change_Dim_Sq = zeros(1,Nd);
for t = 1:Max_gen
    if SF == 1
        %���ɵĽ�
        % ֪����ǰ�ӵڼ�ά��ʼ��֪������,����n����,֪������������,֪�����䷶Χ
        sol = Generate(temp_local_center, Dim_Continuous, Now_Dim_Slide_A, Step_len_exp, Nd, LB, UB);
        fitness = feval(fhd,sol',varargin{:});
        [fitness_min,f_index] = min(fitness);
        if fitness_min <= Fit_min
            Gbest = sol(f_index,:);
            Fit_min = fitness_min;
            temp_local_center = sol(f_index,:); 
            nochange_step_flag = 0; % �ҵ��õľ�����Ϊ0
        else
            % ��ǰά�Ҳ������õľͻ���һά
            Now_Dim_Slide_A = Now_Dim_Slide_A +1;
            % ����������û���ҵ����õľͼ�С����
            nochange_step_flag = nochange_step_flag + 1;
            if nochange_step_flag>Nd/2
                Step_len_exp = Step_len_exp * 2;
                nochange_step_flag = 0;
            end
        end
        if Now_Dim_Slide_A == Nd+1,  Now_Dim_Slide_A = 1;  end
    end
    if SF == 2
        % ״̬2
        %������仯��С��ά������
        [temp,Change_Dim_Sq] = sort(abs(local_history(end,:)-local_history(end-1,:)));

        temp_new_center = local_history(end,:);
        sol = Generate_2(temp_new_center, Step_len_exp, Change_Dim_Sq, Now_Slide_2, Dim_Continuous,s_num, Nd, LB, UB);
        fitness = feval(fhd,sol',varargin{:});
        [fitness_min,f_index] = min(fitness);
        if fitness_min < Fit_min 
            Gbest = sol(f_index,:);
            Fit_min = fitness_min;
        else
            if s_num >= Step_len_exp
                Now_Slide_2 = Now_Slide_2 + 1;
                s_num = 1;
            end
            s_num = s_num + 1;
        end
        if Now_Slide_2 == Nd+1
            Now_Slide_2 = 1;  
            %ת��
            local_history = [local_history;Gbest;];
            Step_len_exp = Step_len_exp * 2;
        end
    end
    if (SF==1) &&  t >Nd*3
        if (fit_save(t-Nd-1)-fit_save(t-1))<1.0e-20
            SF=2;
            %��¼״̬1��ɺ�����Ž�
            local_history = [local_history;Gbest;];
            Step_len_exp = Step_Coefficient*8;
        end
    end    
    
    fit_save(t) = Fit_min;
    if SF == 1
        dim_save(t) = Now_Dim_Slide_A;
    else
        dim_save(t) = Now_Slide_2;
    end
    
end
end

function sol = Generate_2(temp_new_center, Step_len_exp, Change_Dim_Sq, Now_Slide_2, Dim_Continuous,s_num, Nd, LB, UB)
neibor_number = 1;
sol = zeros(Dim_Continuous*2,Nd);
d_change_2 = zeros(1,Dim_Continuous);
for d_c = 1:Dim_Continuous
    temp = Now_Slide_2+d_c-1;
    if temp>Nd
        temp = temp-Nd;
    end
    d_change_2(d_c) = Change_Dim_Sq(temp);
end

for t = 1:Dim_Continuous
    temp_A1 = temp_new_center;
    temp_A2 = temp_new_center;
    d = d_change_2(t);
    base_step1 = (UB(d) - temp_new_center(d))/Step_len_exp;
    base_step2 = (temp_new_center(d) - LB(d))/Step_len_exp;
    temp_A1(d) = temp_A1(d) + base_step1*s_num;
    temp_A2(d) = temp_A2(d) - base_step2*s_num;

    sol(neibor_number,:) = temp_A1;
    neibor_number = neibor_number +1;
    sol(neibor_number,:) = temp_A2;
    neibor_number = neibor_number +1;
end

end
function sol = Generate(A, Dim_Continuous, Now_Dim_Slide_A, Step_len_exp, Nd, LB, UB)
neibor_number = 1;
sol = zeros(Dim_Continuous*2,Nd);
d_change = zeros(1,Dim_Continuous);
for d_c = 1:Dim_Continuous
    d_change(d_c) = Now_Dim_Slide_A+d_c-1;
    if d_change(d_c) >Nd
        d_change(d_c) = d_change(d_c)-Nd;
    end
end
for t = 1:Dim_Continuous
    temp_A1 = A;
    temp_A2 = A;
    d = d_change(t);
    temp_A1(d) = temp_A1(d) + (UB(d)-LB(d))/Step_len_exp;
    temp_A2(d) = temp_A2(d) - (UB(d)-LB(d))/Step_len_exp;
    sol(neibor_number,:) = temp_A1;
    neibor_number = neibor_number +1;
    sol(neibor_number,:) = temp_A2;
    neibor_number = neibor_number +1;
end

end

