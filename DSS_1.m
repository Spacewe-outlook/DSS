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
% 连续几维
Dim_Continuous = Dim_Continuous; %每次迭代计算两个值
% 记录当前移动到第几维
Now_Dim_Slide_A = 1;
% 记录状态2移动到第几维
Now_Slide_2 = 1;
% 计算固定初始位置
center = (LB+UB)/2;
A = (LB+center)/2;
% 保存适应度值
fit_save=zeros(1,Max_gen);
%保存维度变化
dim_save =zeros(1,Max_gen);
% 保存局部最优
local_history = [];
local_history = [local_history;A;];
% 全局最优解
Fit_min = Inf;
% 状态标志
SF = 1;
temp_local_center = A;
s_num = 1;
%维度变化差异排列
Change_Dim_Sq = zeros(1,Nd);
for t = 1:Max_gen
    if SF == 1
        %生成的解
        % 知道当前从第几维开始，知道步长,生成n个解,知道解虚拟中心,知道区间范围
        sol = Generate(temp_local_center, Dim_Continuous, Now_Dim_Slide_A, Step_len_exp, Nd, LB, UB);
        fitness = feval(fhd,sol',varargin{:});
        [fitness_min,f_index] = min(fitness);
        if fitness_min <= Fit_min
            Gbest = sol(f_index,:);
            Fit_min = fitness_min;
            temp_local_center = sol(f_index,:); 
            nochange_step_flag = 0; % 找到好的就重置为0
        else
            % 当前维找不到更好的就换下一维
            Now_Dim_Slide_A = Now_Dim_Slide_A +1;
            % 如果连续多次没有找到更好的就减小步长
            nochange_step_flag = nochange_step_flag + 1;
            if nochange_step_flag>Nd/2
                Step_len_exp = Step_len_exp * 2;
                nochange_step_flag = 0;
            end
        end
        if Now_Dim_Slide_A == Nd+1,  Now_Dim_Slide_A = 1;  end
    end
    if SF == 2
        % 状态2
        %计算出变化最小的维度序列
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
            %转移
            local_history = [local_history;Gbest;];
            Step_len_exp = Step_len_exp * 2;
        end
    end
    if (SF==1) &&  t >Nd*3
        if (fit_save(t-Nd-1)-fit_save(t-1))<1.0e-20
            SF=2;
            %记录状态1完成后的最优解
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

