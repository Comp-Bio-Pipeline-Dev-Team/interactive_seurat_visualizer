
library(ggplot2)
print(paste("ggplot2 version:", packageVersion("ggplot2")))
print(paste("check_linewidth exists:", exists("check_linewidth", where = asNamespace("ggplot2"))))

# Try to load ggtree if present
tryCatch({
  library(ggtree)
  print("ggtree loaded successfully")
}, error = function(e) {
  print(paste("ggtree load failed:", e$message))
})
