shuffle <- function(y, cluster = NULL){
  if(is.null(cluster)){
    cluster <- 1:length(y)
  }
  j <- sample(unique(cluster), length(unique(cluster)), replace = TRUE)
  deck <- 1:length(y)
  deck <- deck[cluster %in% j]
  k <- sample(deck, length(y), replace = TRUE)
  list(deck = deck, k = k)
}

y <- runif(100, 1, 10)
cluster <- sample(c('a','b','c','d','e','f','g'), 100, replace = TRUE)

a <- shuffle(y, cluster)
data.frame(y[a$k], cluster[a$k])
