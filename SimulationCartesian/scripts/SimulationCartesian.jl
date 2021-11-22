using DataFrames
using Distributions
import Random
using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
import RDatasets

#Fixed parameters
MOQ=100
VF=0.5
LTF=0.5
MOC=1 #given by the problem statement

Nweeks=8
mult=1 #multiplier to give more importance to the missing orders

#Leadtime
leadTime=[5,5,3,3,4,6,5,2,5,3,4,4,3,4,5,4,2,2,4,4,3,5,4,3,4,4,4,5,3,4,5,4,5,5,3,5,5,4,2]
lowerLT=minimum(leadTime)
upperLT=maximum(leadTime)
pdLeadTime=TriangularDist(lowerLT, upperLT)

#Creation of a distribution function to make a random demand
sales=[108 80 104 100 112 104 100 80 108 92 112 92 112 92 80 96 80 112 96 88 104 88 104 88 80 100 96 92 108 84 80 112];
lower=minimum(sales)
upper=maximum(sales)
pd=TriangularDist(lower, upper)

    

function evaluate(ind::CGPInd)
    accuracy = 0.0
    #input=[LT, ADU, MOQ, VF,LTF,MOC,0,1]
    input=[rand(pdLeadTime),mean(rand(pd)+rand(pd)+rand(pd)+rand(pd)+rand(pd)),MOQ,VF,LTF,MOC,0,1]
    surplus = 0
    missed = 0

    #Normalization
    normalized=[input[1]/(2*upperLT), input[2]/(4*upper),input[3],input[4],input[5],input[6]]

    stock0=80

    for i in 1:Nweeks
        demand=[]
        for i in 1:5
            push!(demand,round(Int64, rand(pd)))
        end 
        
        #Calculation of the buffer levels
        out = process(ind, normalized)

        # Buffer level calculation 
        out[1] =round(Int64, (out[1] + 1) * MOQ) 
        out[2] = round(Int64, (out[2] + 1) * MOQ + out[1])  
        out[3] =round(Int64,(out[3] + 1) * MOQ + out[2]) 

        #Simulation of the week
        for i in 1:5
            if demand[i]>stock0
                order=out[3]
                stock=order-(demand[i]-stock0)
                if stock<out[2]
                    order=out[3]+(out[3]-stock)
                    stock=out[3]
                end
            elseif demand[i]<stock0
                    order=out[3]
                    stock=out[3]
            end

            #fitness Calculation
            if stock >0
                accuracy=accuracy+out[3]-stock # surplus
                surplus = surplus + out[3] - stock
            elseif stock<0
                accuracy=accuracy+stock*mult # missed orders
                missed = missed + stock
                stock=0
            end
            stock0=stock
        end
 
        #Calculation of the new input parameters
        input=[round(Int64, rand(pdLeadTime)),mean(demand),0,VF,LTF,MOC,0,1]

        #Normalization
        normalized=[input[1]/(2*upperLT), input[2]/(4*upper),input[3],input[4],input[5],input[6]]

     end
     [accuracy,-surplus,missed]
end

cfg = get_config("cfg/SimulationCartesian.yaml")
fit(i::CGPInd) = evaluate(i) #compute the fitnes of an individual with the function evaluate
i = CGPInd(cfg)#create  a new individual
inputs = rand(cfg.n_in)
println("Inputs: ", inputs)
out = process(i, inputs)
println("Outputs: ", out)
fitness = evaluate(i)
println(fitness)

#mutate(i::CGPInd) = goldman_mutate(cfg, i)
#e = CGPEvolution(cfg, fit)
#run!(e)

