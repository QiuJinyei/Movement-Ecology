library(moveHMM)
library(bayesmove)
library(readxl)

x <-read_xls("C:/Users/Administrator/Desktop/wash/23w/23w_f.xls")
prep_data <- prepData(x,type="LL",coordNames=c("location_long","location_lat"))
write.csv(prep_data,"C:/Users/Administrator/Desktop/wash/23w/23w_m.csv")
head(prep_data)
# 绘制步长的直方图
breaks <- c(0, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1, 2,3)
hist(prep_data$step,breaks = breaks)
# 步长参数的起始值
stepMean0 <- c(0.03,0.1653) # initial means (one for each state)
stepSD0 <- c(0.027,0.1452) # initial standard deviations (one for each state)
zeromass0 <- c(0.003,0) # step zero-mass
stepPar0 <- c(stepMean0, stepSD0,zeromass0)
# 长度为零的步的索引值
whichzero <- which(prep_data$step == 0)
# 数据集中长度为零的步长所占的比例
length(whichzero)/nrow(prep_data)
## [1] 0
# 绘制转向角的直方图
hist(prep_data$angle, breaks = seq(-pi, pi, length = 15))
# 转向角参数的起始值
angleMean0 <- c(3.107, -0.053) # initial means (one for each state)
angleCon0 <- c(0.435,0.191) # initial concentrations (one for each state)
anglePar0 <- c(angleMean0, angleCon0)
#检查哪些行包含 NA 值
rows_with_na <- which(rowSums(is.na(prep_data)) > 0)
# 删除含有 NA 的行
prep_data_clean <- prep_data[-rows_with_na, ]
# 现在使用 prep_data_complete 来进行模型拟合
m <- fitHMM(data = prep_data, nbStates = 2, stepPar0 = stepPar0, anglePar0 = anglePar0,stepDist="gamma",angleDist="vm"
            ,formula = ~1)
m
AIC(m)
# For reproducibility
set.seed(12345)
# Number of tries with different starting values
niter <- 25
allm <- list()


for(i in 1:niter) {
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(0.01, 0.1),
                     max = c(0.1, 0.3))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(0.01, 0.1),
                   max = c(0.1, 0.3))
  zeromass0 <- c(0.01,0.001)
  # Turning angle mean
  angleMean0 <- c(pi,0)
  # Turning angle concentration
  angleCon0 <- runif(2,
                     min = c(0.05, 0.5),
                     max = c(0.2, 2))
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0,zeromass0)
  anglePar0 <- c(angleMean0, angleCon0)
  allm[[i]] <- fitHMM(data = prep_data, nbStates = 2, stepPar0 = stepPar0,
                      anglePar0 = anglePar0)
}
# Extract likelihoods of fitted models
allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
allnllk
#[1] 581.8873 581.8873 581.8873 581.8873 581.8873 581.8873 581.8873 581.8873 602.8307 581.8873 581.8873
#[12] 602.8307 581.8873 581.8873 581.8873 581.8873 581.8873 602.8307 581.8873 581.8873 581.8873 581.8873
#[23] 581.8873 581.8873 581.8873
# Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)
# Best fitting model
mbest <- allm[[whichbest]]
mbest

AIC(mbest)
CI(mbest)
plot(mbest, plotCI=TRUE)
states <- viterbi(mbest)
states<-data.frame(states)
prep_data$states = states
sp <- stateProbs(mbest)
head(sp)
plotStates(mbest)
# compute the pseudo-residuals
pr <- pseudoRes(mbest)
pr<-data.frame(pr)
prep_data$pr = pr
# time series, qq-plots, and ACF of the pseudo-residuals
plotPR(mbest)
#Computing pseudo-residuals... Note: Some angles are equal to pi, and the corresponding pseudo-residuals are not included
#DONE
#Error in plot.new() : figure margins too large


#______________________________________________________________________________________
# 加载所需包
library(terra)
library(raster)

# 读取 DEM 地图
dem <- rast("C:/Users/Administrator/Desktop/TPG/TPGGC/tpgdem.tif")

# 创建包含 x, y 和 states 列的数据框
lldata <- data.frame(x = prep_data$x, y = prep_data$y, states = prep_data$states)

# 将点坐标转换为空间对象，并包含 states 列
points <- vect(lldata, geom = c("x", "y"), crs = crs(dem))

# 获取点的范围
points_extent <- ext(points)

# 扩大范围以便更好地显示周围区域（例如，扩大10%）
buffer <- 0.1 * (points_extent[2] - points_extent[1])
new_extent <- ext(points_extent[1] - buffer, points_extent[2] + buffer,
                  points_extent[3] - buffer, points_extent[4] + buffer)

# 裁剪 DEM 地图
dem_cropped <- crop(dem, new_extent)

# 调整图形参数，增加右边距
par(mar = c(5, 4, 4, 10) + 0.1)


# 创建自定义白色到黄色的颜色梯度
custom_colors <- colorRampPalette(c("#ffffff", "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c"))(10)

# 绘制裁剪后的 DEM 地图，使用自定义白色到黄色的颜色梯度渲染高程值，并去掉颜色条
plot(dem_cropped, col = custom_colors, main = "Male Behavioral States", legend = FALSE)

# 获取所有唯一的状态
unique_states <- unique(points$states)

# 为每个状态分配颜色
state_colors <- rainbow(length(unique_states))
# 确保有两个状态，如果状态超过两个，请根据需要调整颜色分配
state_colors <- c("red", "darkblue")

# 创建一个状态到颜色的映射
state_color_map <- setNames(state_colors, unique_states)

# 根据 states 列的值设置颜色
point_colors <- state_color_map[as.character(points$states)]

# 添加动物轨迹点，颜色根据 states 改变
plot(points, add = TRUE, col = point_colors, pch = 1, cex = 0.5)

# 添加图例，并设置黑色边框
legend("topright", inset = c(-0.1, 0), legend = unique_states, col = state_colors, pch = 1, title = "States",
       box.lwd = 2, box.col = "black", xpd = TRUE)

#------------------------------------------------------------------------------------

# 加载所需包
library(terra)
library(raster)

# 读取底图
base_map <- rast("C:/Users/Administrator/Desktop/TPG/TPGSATE/TPG14.tif")

# 创建包含 x, y 和 states 列的数据框
lldata <- data.frame(x = prep_data$x, y = prep_data$y, states = prep_data$states)

# 将点坐标转换为空间对象，并包含 states 列，假设原始坐标系是 WGS84
points <- vect(lldata, geom = c("x", "y"), crs = "EPSG:4326")

# 将点数据转换到底图的坐标系
points <- project(points, crs(base_map))

# 获取点的范围
points_extent <- ext(points)

# 打印点的范围和底图的范围以进行检查
print(points_extent)
print(ext(base_map))

# 检查点的范围是否与底图的范围重叠
intersection_extent <- intersect(points_extent, ext(base_map))
if (is.null(intersection_extent)) {
  stop("点的范围与底图的范围不重叠，请检查数据。")
}

# 扩大范围以便更好地显示周围区域（例如，扩大10%）
buffer <- 0.1 * (points_extent[2] - points_extent[1])
new_extent <- ext(points_extent[1] - buffer, points_extent[2] + buffer,
                  points_extent[3] - buffer, points_extent[4] + buffer)

# 打印新的裁剪范围以进行检查
print(new_extent)

# 确保裁剪范围在底图范围内
new_extent <- intersect(new_extent, ext(base_map))

# 打印新的裁剪范围以进行检查
print(new_extent)

# 确保 new_extent 是 SpatExtent 对象
if (!inherits(new_extent, "SpatExtent")) {
  stop("new_extent 不是 SpatExtent 对象，请检查数据。")
}

# 裁剪底图
base_map_cropped <- crop(base_map, new_extent)

# 绘制裁剪后的底图
plot(base_map_cropped, main = "Behavior States")

# 获取所有唯一的状态
unique_states <- unique(points$states)

# 为每个状态分配颜色
state_colors <- rainbow(length(unique_states))

# 创建一个状态到颜色的映射
state_color_map <- setNames(state_colors, unique_states)

# 根据 states 列的值设置颜色
point_colors <- state_color_map[as.character(points$states)]

# 添加动物轨迹点，颜色根据 states 改变
plot(points, add = TRUE, col = point_colors, pch = 1, cex = 0.5)

# 添加图例
legend("topright", legend = unique_states, col = state_colors, pch = 1, title = "States")
#See ?register_google for details.
write.csv(states,file="C:/Users/Administrator/Desktop/wash/23w/states23wf",row.names=TRUE,quote=TRUE)
write.csv(prep_data,file="C:/Users/Administrator/Desktop/wash/23w/states23wf.csv",row.names = TRUE,quote = TRUE)
