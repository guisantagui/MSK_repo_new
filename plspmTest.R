
if(!require(tester)) install.packages("tester")
if(!require(turner)) install.packages("turner")
if(!require(diagram)) install.packages("diagram")
if(!require(shape)) install.packages("shape")
if(!require(amap)) install.packages("amap")
if(!require(plspm)) install.packages(pkgs = "C:/Users/Guillem/Documents/R/win-library/plspm_0.4.9.tar.gz", repos = NULL, type = "source")
library(plspm)

data("spainfoot")

head(spainfoot)

# rows of the inner model matrix
Attack = c(0, 0, 0) 
Defense = c(0, 0, 0) 
Success = c(1, 1, 0)
# path matrix created by row binding 
foot_path = rbind(Attack, Defense, Success)
# add column names (optional) 
colnames(foot_path) = rownames(foot_path)

innerplot(foot_path)

foot_blocks = list(1:4, 5:8, 9:12)

foot_modes = c("A", "A", "A")

foot_modes2 = c("A", "A", "B")

foot_pls = plspm(spainfoot, foot_path, foot_blocks, modes = foot_modes)

foot_pls
