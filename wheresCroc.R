
myFunction<- function(moveInfo, readings, positions, edges, probs) {
  
  mem <- moveInfo$mem
  
  W <- 1
  n <- 40
  
  # ---------------------------------
  # HELPER: SAFE NORMALISATION
  # ---------------------------------
  
  normalise <- function(x, fallback) {
    
    sx <- sum(x)
    
    if (is.finite(sx) && sx > 0) {
      return(x / sx)
    }
    
    sf <- sum(fallback)
    
    if (is.finite(sf) && sf > 0) {
      return(fallback / sf)
    }
    
    # Last-resort fallback if no valid distribution exists
    rep(1 / length(x), length(x))
  }
  
  # ---------------------------------
  # 1. INITIALISE NEW GAME
  # ---------------------------------
  
  if (mem$status < 2) {
    
    # Static setup: only performed once
    if (mem$status == 0) {
      
      # Adjacency list
      adj <- vector("list", n)
      
      for (r in seq_len(nrow(edges))) {
        
        a <- edges[r, 1]
        b <- edges[r, 2]
        
        adj[[a]] <- c(adj[[a]], b)
        adj[[b]] <- c(adj[[b]], a)
      }
      
      adj <- lapply(adj, function(x) sort(unique(x)))
      
      # Transition matrix
      # Croc can stay or move to any neighbouring waterhole
      trans <- matrix(0, n, n)
      
      for (i in 1:n) {
        
        options <- c(i, adj[[i]])
        
        trans[i, options] <- 1 / length(options)
      }
      
      # BFS distances between all waterholes
      dist <- matrix(Inf, n, n)
      
      for (s in 1:n) {
        
        dist[s, s] <- 0
        
        queue <- s
        head <- 1
        
        while (head <= length(queue)) {
          
          v <- queue[head]
          head <- head + 1
          
          for (w in adj[[v]]) {
            
            if (is.infinite(dist[s, w])) {
              
              dist[s, w] <- dist[s, v] + 1
              
              queue <- c(queue, w)
            }
          }
        }
      }
      
      # Precompute candidate two-move plans
      plans <- vector("list", n)
      
      for (p in 1:n) {
        
        opts1 <- c(0, p, adj[[p]])
        
        out <- list()
        
        for (m1 in opts1) {
          
          pos1 <- if (m1 == 0) p else m1
          
          opts2 <- c(0, pos1, adj[[pos1]])
          
          for (m2 in opts2) {
            
            pos2 <- if (m2 == 0) pos1 else m2
            
            searched <- c(
              if (m1 == 0) p,
              if (m2 == 0) pos1
            )
            
            out[[length(out) + 1]] <- list(
              moves = c(m1, m2),
              searched = unique(searched),
              end = pos2
            )
          }
        }
        
        plans[[p]] <- out
      }
      
      # Store static information
      mem$adj <- adj
      mem$trans <- trans
      mem$dist <- dist
      mem$plans <- plans
      
      # Possible Croc positions after one movement
      mem$reach <- lapply(
        1:n,
        function(p) c(p, adj[[p]])
      )
    }
    
    # Per-game initialisation
    mem$state <- rep(1 / n, n)
    mem$first <- TRUE
    mem$status <- 2
  }
  
  # Retrieve stored information
  adj <- mem$adj
  trans <- mem$trans
  dist <- mem$dist
  plans <- mem$plans
  reach <- mem$reach
  
  state <- mem$state
  
  # ---------------------------------
  # 2. PREDICT CROCODILE LOCATION
  # ---------------------------------
  
  if (!mem$first) {
    
    state <- as.vector(state %*% trans)
  }
  
  # ---------------------------------
  # 3. SENSOR UPDATE
  # ---------------------------------
  
  # Calculate log-likelihoods to avoid numerical underflow
  log_lik <-
    
    dnorm(
      readings[1],
      probs$salinity[, 1],
      probs$salinity[, 2],
      log = TRUE
    ) +
    
    dnorm(
      readings[2],
      probs$phosphate[, 1],
      probs$phosphate[, 2],
      log = TRUE
    ) +
    
    dnorm(
      readings[3],
      probs$nitrogen[, 1],
      probs$nitrogen[, 2],
      log = TRUE
    )
  
  # Convert log-likelihoods into relative likelihoods
  max_log <- max(log_lik)
  
  if (is.finite(max_log)) {
    
    lik <- exp(log_lik - max_log)
    
  } else {
    
    # If all likelihoods are impossible, ignore this sensor reading
    lik <- rep(1, n)
  }
  
  prior <- state
  
  state <- normalise(state * lik, prior)
  
  # ---------------------------------
  # 4. BACKPACKER UPDATE
  # ---------------------------------
  
  me <- positions[3]
  
  before <- state
  
  eaten_at <- NA
  
  valid <- rep(TRUE, n)
  
  for (k in 1:2) {
    
    p <- positions[k]
    
    if (is.na(p)) next
    
    if (p < 0) {
      
      eaten_at <- -p
      
    } else {
      
      # A living backpacker confirms this location is safe
      valid[p] <- FALSE
    }
  }
  
  # Initially Croc cannot be at the ranger's starting position
  if (mem$first) {
    
    valid[me] <- FALSE
  }
  
  if (!is.na(eaten_at)) {
    
    # A backpacker being eaten identifies Croc's location
    state <- rep(0, n)
    
    state[eaten_at] <- 1
    
  } else {
    
    # Remove impossible locations
    state[!valid] <- 0
    
    # Keep fallback probabilities consistent with the evidence
    fallback <- before
    
    fallback[!valid] <- 0
    
    if (sum(fallback) <= 0) {
      
      fallback <- as.numeric(valid)
    }
    
    state <- normalise(state, fallback)
  }
  
  mem$first <- FALSE
  
  # ---------------------------------
  # 5. CHOOSE RANGER MOVES
  # ---------------------------------
  
  target <- which.max(state)
  
  best_score <- -Inf
  best <- NULL
  
  near <- dist[me, target] <= 1
  
  for (pl in plans[[me]]) {
    
    if (near) {
      
      # Search the most likely location this turn
      if (!(target %in% pl$searched)) next
      
      # Remove locations searched during this plan
      rest <- state
      
      rest[pl$searched] <- 0
      
      # Predict Croc's next position
      nxt <- as.vector(rest %*% trans)
      
      # Evaluate the best possible follow-up search
      score <-
        
        max(nxt[reach[[pl$end]]]) * W +
        
        sum(state[pl$searched])
      
    } else {
      
      # Move closer to the most likely Croc location
      score <-
        
        -dist[pl$end, target] -
        
        1e-3 * length(pl$searched)
    }
    
    if (score > best_score) {
      
      best_score <- score
      
      best <- pl
    }
  }
  
  move1 <- best$moves[1]
  move2 <- best$moves[2]
  
  # ---------------------------------
  # 6. UPDATE BELIEF AFTER SEARCH
  # ---------------------------------
  
  # Croc cannot be at locations searched unsuccessfully
  state[best$searched] <- 0
  
  # Construct a fallback distribution over unsearched locations
  remaining <- rep(TRUE, n)
  
  remaining[best$searched] <- FALSE
  
  state <- normalise(
    state,
    as.numeric(remaining)
  )
  
  # ---------------------------------
  # 7. RETURN
  # ---------------------------------
  
  mem$state <- state
  
  moveInfo$mem <- mem
  
  moveInfo$moves <- c(move1, move2)
  
  return(moveInfo)
}