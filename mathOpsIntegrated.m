function [added, multed] = mathOpsIntegrated(in1, in2)
added = 0;

coder.cinclude('adder.h');
added = coder.ceval('adder', in1, in2);
multed = in1*in2;
end