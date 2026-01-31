
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Set target version
BiocManager::install(version = "3.19", ask = FALSE, update = FALSE)

max_attempts <- 10
attempt <- 1

while (attempt <= max_attempts) {
    print(paste("Attempt", attempt, "to clean 'too new' packages..."))
    
    invalid <- BiocManager::valid(fix=FALSE)
    
    if (isTRUE(invalid)) {
        print("Environment is valid!")
        break
    }
    
    if (is.null(invalid$too_new)) {
        print("No 'too new' packages found. Environment is clean enough.")
        break
    }
    
    pkgs_to_remove <- rownames(invalid$too_new)
    print(paste("Removing", length(pkgs_to_remove), "packages..."))
    remove.packages(pkgs_to_remove)
    
    attempt <- attempt + 1
}

if (attempt > max_attempts) {
    print("Warning: Max attempts reached. Continuing anyway...")
}

print("Installing clusterProfiler (version 3.19)...")
BiocManager::install("clusterProfiler", version="3.19", ask=FALSE, force=TRUE, update=TRUE)
