clear all
close all
load("proj_fit_35.mat")

X1id = id.X{1};
X2id = id.X{2};
Yid = id.Y;

X1val = val.X{1};
X2val = val.X{2};
Yval = val.Y;

n_id = id.dims(1);
n_val = val.dims(1);

m = 30; %grad configurabil

vecMSE_id = zeros(1, m);
vecMSE_val = zeros(1, m);

for g = 1 : m
    nr = ((g+2)*(g+1))/2;
    phi_id = zeros(n_id*n_id, nr);
    phi_val = zeros(n_val*n_val, nr);
    k = 1;
    for i = 1 : n_id 
        for j = 1 : n_id
            p = 0;
            q = 0;
            for n = 1 : nr
                phi_id(k, n) = (X1id(i).^ (p - q)).* (X2id(j).^ q);
                if(q == p && p < g)
                    p = p + 1;
                    q = 0;
                else
                    q = q + 1;
                end         
             end
             k = k + 1;
        end   
    end
    
    k = 1;
    Y_id = zeros(n_id*n_id, 1);
    for i = 1 : n_id
        for j = 1 : n_id
            Y_id(k, 1) = Yid(i, j);
            k = k + 1;
        end
    end

    theta = phi_id \ Y_id;

    Yid_hat = phi_id * theta;
    
    YY_id = zeros (n_id, n_id);
    k = 1;
    for i = 1 : n_id
        for j = 1 : n_id
            YY_id(i, j) = Yid_hat(k, 1);
            k = k + 1;
        end
    end
    
    s_id = 0;
    for i = 1 : n_id
        s_id = s_id + (Yid(i) - YY_id(i))^2;
    end
    MSE_id = 1/n_id * s_id;
    vecMSE_id(g) = MSE_id;
    
    k = 1;
    for i = 1 : n_val 
        for j = 1 : n_val
            p = 0;
            q = 0;
            for n = 1 : nr
                phi_val(k, n) = (X1val(i).^ (p - q)).* (X2val(j).^ q);
                if(q == p && p < g)
                    p = p + 1;
                    q = 0;
                else
                    q = q + 1;
                end         
             end
             k = k + 1;
        end   
    end
    
    k = 1;
    Y_val = zeros(n_val*n_val, 1);
    for i = 1 : n_val
        for j = 1 : n_val
            Y_val(k, 1) = Yval(i, j);
            k = k + 1;
        end
    end
  
    Yval_hat = phi_val * theta;
    
    k = 1;
    YY_val = zeros (n_val, n_val);
    for i = 1 : n_val
        for j = 1 : n_val
            YY_val(i, j) = Yval_hat(k, 1);
            k = k + 1;
        end
    end
               
    s_val = 0;
    for i = 1 : n_val
        s_val = s_val + (Yval(i) - YY_val(i))^2;
    end
    MSE_val = 1/n_val * s_val;
    vecMSE_val(g) = MSE_val;

    if(g == 5)
        subplot(2, 1, 1)
        mesh(X1val, X2val, Yval)
        title('Date validare')
        xlabel('X1')
        ylabel('X2')
        zlabel('Y')
        
        subplot(2, 1, 2)
        mesh(X1val, X2val, YY_val)
        title('Aproximare date validare MSE minim')
        xlabel('X1')
        ylabel('X2')
        zlabel('Y')
        figure
    end
end

subplot(2, 1, 1)
mesh(X1id, X2id, Yid)
title('Date identificare')
xlabel('X1')
ylabel('X2')
zlabel('Y')

subplot(2, 1, 2)
mesh(X1val, X2val, YY_val)
title('Aproximare date validare')
xlabel('X1')
ylabel('X2')
zlabel('Y')
figure

plot(1 : m, vecMSE_id)
title('Valori MSE identificare')
xlabel('Gradul polinomului m')
ylabel('MSE')
figure

plot(1 : m, vecMSE_val)
title('Valori MSE validare')
xlabel('Gradul polinomului m')
ylabel('MSE')

minim_MSEval = 100000;
for i = 1 : m
    if vecMSE_val(i) < minim_MSEval
        minim_MSEval = vecMSE_val(i);
    end
end