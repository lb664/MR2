# 135 characters #####################################################################################################################
#' @title Sampler of decomposable graphs
#' @description Internal functions that perform one Markov chain Monte Carlo iteration to sample from the posterior distribution of a 
#' decomposable graph. The R code of this function is a modified version of the Matlab code in \insertCite{Talhouk2012;textual}{MR2}. 
#' For details, see \insertCite{Alexopoulos2021;textual}{MR2}
#'
#' @references
#' \insertAllCited{}


########## Sample_G ##########

Sample_G <- function(W, GCurr, samplerPar, delta = 2)
{
  
  burnin <- samplerPar$GPar$burnin
  niter <- samplerPar$GPar$niter
  
  n <- dim(W)[1]
  J <- dim(W)[2]
  Phi <- diag(rep(1, J))
  WW <- t(W) %*% W
  
  G <- generate_graph_zo(GCurr, delta, Phi, n, WW, burnin, niter)[[1]]
  output <- G
  
  return(output)
}

########## ripcliques_to_jtree_zo ##########

ripcliques_to_jtree_zo <- function(cliques)
{
  p <- dim(cliques)[1]
  clique_sizes <- rep(0, p)
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  score <- rep(0, p)
  jtree <- matrix(0, p, p)
  if (num_cliques >= 2)
  {
    for (i in 2 : num_cliques)
    {
      for (k in 1 : (i - 1))
      {
        score[k] <- sum(intersect_zo(cliques[, i], cliques[, k]))
      }
      if (max(score) != 0)
      {
        argmax <- which(score == max(score))[1]
        j <- argmax
        jtree[i, j] <- 1
        jtree[j, i] <- 1
      }
    }
  }
  output <- jtree
  
  return(output)
}

########## sepsize_zo ##########

sepsize_zo <- function(cliques, jtree)
{
  p <- dim(cliques)[1]
  sepsize <- matrix(0, p, p)
  clique_sizes <- rep(0, p)
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  if (num_cliques >= 1)
  {
    for (i in 1 : num_cliques)
    {
      if ((i + 1) <= num_cliques)
      {
        for (j in (i + 1) : num_cliques)
        {
          if (jtree[i,j] == 1)
          {
            sepsize[i, j] <- sum(intersect_zo(cliques[, i], cliques[, j]))
            sepsize[j, i] <- sepsize[i, j]
          }
        }
      }
    }
  }
  output <- sepsize
  
  return(output)
}

########## reachability_graph ##########

reachability_graph <- function(g)
{
  p <- dim(g)[1]
  A <- g
  reach_graph <- matrix(0, p, p)
  for (i in 1 : (p - 1))
  {
    reach_graph <- reach_graph + A
    A <- A %*% g
  }
  reach_graph[which(reach_graph > 0)] <- 1
  output <- reach_graph
  
  return(output)
}

########## next_edge_candidate ##########

next_edge_candidate <- function(g_current)
{
  p <- max(dim(g_current)[1], dim(g_current)[2])
  iota <- sample(p)
  i <- iota[1]
  j <- iota[2]
  edge_candidate <- c(i, j)
  output <- edge_candidate
  
  return(output)
}

########## check_edge_delete_zo ##########

check_edge_delete_zo <- function(i, j, cliques)
{
  p <- dim(cliques)[1]
  delete_ok <- 1
  C <- rep(0, p)
  count_cliques <- 0
  clique_sizes <- rep(0, p)
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  edge_col_vec <- rep(0, p)
  edge_col_vec[i] <- 1
  edge_col_vec[j] <- 1
  if (num_cliques >= 1)
  {
    for (k in 1 : num_cliques)
    {
      clique_k <- cliques[, k]
      if (sum(intersect_zo(edge_col_vec, clique_k)) == 2)
      {
        C <- clique_k
        count_cliques <- count_cliques + 1
        if (count_cliques > 1) 
        {
          delete_ok <- 0
          C <- rep(9, 3)
          
          break
        }
      }
    }
  }
  output <- list(delete_ok, C)
  
  return(output)
}

########## check_edge_delete_zo ##########

find_clique_containing_zo<-function(a, cliques){

  p <- dim(cliques)[1]
  all_indicies <- rep(0, p)
  all_indicies <- which(cliques[a, ]!=0)
  index_a <- min(all_indicies)
  output <- index_a
  
  return(output)
}

########## neighbours_vertex_zo ##########

neighbours_vertex_zo <- function(g, v_i)
{
  p <- dim(g)[1]
  ns <- rep(0, p)
  ns <- g[, v_i]
  output <- ns
  
  return(output)
}

########## parents_vertex_zo ##########

parents_vertex_zo <- function(g, i)
{
  p <- dim(g)[1]
  ps <- rep(0,p)
  nbs <- neighbours_vertex_zo(g, i);
  ps <- nbs
  ps[i : p] <- 0
  output <- ps
  
  return(output)
}

########## check_edge_add_same_component_zo ##########

check_edge_add_same_component_zo <- function(a, b, jtree, sepsize, cliques)
{
  edge_ok <- 0
  quit <- 0
  fork <- 0
  CASE_same_branch <- 0
  CASE_sats_on_b2_branch <- 0
  locate_a2 <- 0
  p <- dim(cliques)[1]
  ab_edge_vec <- rep(0, p)
  ab_edge_vec[a] <- 1
  ab_edge_vec[b] <- 1
  index_b <- find_clique_containing_zo(b, cliques)
  index_a <- find_clique_containing_zo(a, cliques)
  index_b2 <- max(c(index_a, index_b))
  if (index_b2 == index_b)
  {
    b2 <- b
    a2 <- a
  }else{
    b2 <- a
    a2 <- b
  }
  index <- index_b2
  next_index <- which(parents_vertex_zo(jtree, index) !=0 )
  if ((index_b2 - 1) >= 1)
  {
    for (dummy in 1 : (index_b2 - 1))
    {
      if (sum(next_index) == 0 & locate_a2 == 0)
      {
        CASE_same_branch <- 0
        
        break
      }
      clique_next_index <- (cliques[, next_index])
      if (clique_next_index[a2] == 1 & sum(next_index) != 0)
      {
        index_a2 <- next_index
        locate_a2 <- 1
        CASE_same_branch <- 1
        
        break
      }
      index <- next_index
      next_index <- which(parents_vertex_zo(jtree, index) != 0)
    }
  }
  if (CASE_same_branch == 0)
  {
    index_a2 <- find_clique_containing_zo(a2, cliques)
  }
  clique_a2_int_clique_b2 <- intersect_zo(cliques[, index_a2], cliques[, index_b2])
  s <- sum(clique_a2_int_clique_b2)
  if (s == 0)
  {
    edge_ok <- 0
    C <- rep(0, p)
    quit <- 1
  }
  
  if (jtree[index_a2, index_b2] == 1)
  {
    edge_ok <- 1
    C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
    quit <- 1
  }
  if (quit == 0)
  {
    if (CASE_same_branch == 1)
    {
      bottom <- index_b2
      next_parent_b2 <- which(parents_vertex_zo(jtree, index_b2) != 0)
      while ((sum(next_parent_b2) > 0) & (next_parent_b2 >= index_a2))
      {
        if (sepsize[next_parent_b2, bottom] == s)
        {
          edge_ok <- 1
          C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
          quit <- 1
          
          break
        }
        bottom <- next_parent_b2
        next_parent_b2 <- which(parents_vertex_zo(jtree, bottom) !=0 )
        if (length(next_parent_b2) == 0)
        {
          
          break
        }
      }
    }else{
      next_parent_a2 <- which(parents_vertex_zo(jtree, index_a2) != 0)
      next_parent_b2 <- which(parents_vertex_zo(jtree, index_b2) != 0)
      ancestors_a2 <- rep(0, index_a2)
      ancestors_b2 <- rep(0, index_b2)
      ancestors_a2_col_vec <- rep(0, p)
      ancestors_b2_col_vec <- rep(0, p)
      next_ancestor_a2 <- index_a2
      next_ancestor_b2 <- index_b2
      if (index_a2 >= 1)
      {
        for (count in 1 : index_a2)
        {
          ancestors_a2[count] <- next_ancestor_a2
          ancestors_a2_col_vec[next_ancestor_a2] <- 1
          next_ancestor_a2 <- which(parents_vertex_zo(jtree, next_ancestor_a2) != 0)
          if (sum(next_ancestor_a2) == 0)
          {
            
            break
          }
        }
      }
      if (index_b2 >= 1)
      {
        for (count in 1 : index_b2)
        {
          ancestors_b2[count] <- next_ancestor_b2
          ancestors_b2_col_vec[next_ancestor_b2] <- 1
          next_ancestor_b2 <- which(parents_vertex_zo(jtree, next_ancestor_b2) != 0)
          if (sum(next_ancestor_b2) == 0)
          {
            
            break
          }
        }
      }
      fork_set <- intersect_zo(ancestors_a2_col_vec, ancestors_b2_col_vec)
      fork <- max(which((fork_set != 0)))
      num_seps_branch_b2 <- length(which(ancestors_b2 > fork))
      if (num_seps_branch_b2 > 0)
      {
        for (count in 1 : num_seps_branch_b2)
        {
          index_row <- ancestors_b2[count]
          index_col <- ancestors_b2[count + 1]
          if (sepsize[index_row, index_col] == s)
          {
            edge_ok <- 1
            C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
            quit <- 1
            CASE_sats_on_b2_branch <- 1
            
            break
          }
        }
      }
      if (quit == 0 & CASE_sats_on_b2_branch == 0)
      {
        num_seps_branch_a2 <- length(which(ancestors_a2 > fork))
        if (num_seps_branch_a2 > 0)
        {
          for (count in 1 : num_seps_branch_a2)
          {
            index_row <- ancestors_a2[count + 1]
            index_col <- ancestors_a2[count]
            if (sepsize[index_row, index_col] == s)
            {
              edge_ok <- 1
              C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
              quit <- 1
              
              break
            }
          }
        }
      }
    }
  }
  if (quit == 0)
  {
    edge_ok <- 0
    C <- NULL
  }
  output <- list(edge_ok, C)
  
  return(output)
}

########## next_graph_candidate_zo ##########

next_graph_candidate_zo <- function(g_current, jtree, sepsize, cliques, reach_graph)
{
  p <- max(dim(g_current)[1], dim(g_current)[2])
  g_proposal <- g_current
  CASE_add <- 0
  CASE_delete <- 0
  while(identical(g_proposal, g_current))
  {
    edge <- next_edge_candidate(g_current)
    i <- edge[1]
    j <- edge[2]
    if (g_current[i, j] == 1)
    {
      ss <- check_edge_delete_zo(i, j, cliques)
      CASE_delete <- ss[[1]]
      C_potential_delete <- ss[[2]]
      if (CASE_delete == 1)
      {
        g_proposal[i, j] <- 0
        g_proposal[j, i] <- 0
        CASE_delete <- 1
        a <- i
        b <- j
        C <- 
        C_potential_delete
      }
    }else{
      if (reach_graph[i, j] == 0)
      {
        CASE_add <- 1
        a <- i
        b <- j
        C <- c(a, b)
        g_proposal[i, j] <- 1
        g_proposal[j, i] <- 1
      }else{
        ss <- check_edge_add_same_component_zo(i, j, jtree, sepsize, cliques)
        CASE_add <- ss[[1]]
        C_potential_add <- ss[[2]]
        if (CASE_add == 1)
        {
          a <-i
          b <- j
          C <- C_potential_add
          g_proposal[i, j] <- 1
          g_proposal[j, i] <- 1
        }
      }
    }
  }
  output <- list(g_proposal, a, b, C, CASE_add, CASE_delete)
  
  return(output)
}

########## ln_h_ratio_zo ##########

ln_h_ratio_zo <- function(C, a, b, delta, Phi)
{
  C_nodes <- which(C != 0)
  D <- c(a, b)
  Sq2 <- setdiff(C_nodes, D)
  numSq2 <- length(Sq2)
  delta_star <- delta - 1
  Phi_CC <- matrix(0, length(C_nodes), length(C_nodes))
  Phi_CC <- rbind(cbind(matrix(Phi[Sq2, Sq2], length(Sq2), length(Sq2)), matrix(Phi[Sq2,D], length(Sq2), length(D))), 
                  cbind(matrix(Phi[D, Sq2], length(D), length(Sq2)), matrix(Phi[D, D], length(D), length(D))))
  L <- t(chol(Phi_CC))
  L_DD <- L[c(numSq2 + 1, numSq2 + 2), c(numSq2 + 1, numSq2 + 2)]
  l_aa <- L_DD[1, 1]
  l_bb <- L_DD[2, 2]
  l_ab <- L_DD[2, 1]
  ln_top <- (delta_star + numSq2 + 2) * (log(l_aa) + log(l_bb))
  ln_bottom <- (delta_star + numSq2 + 1) / 2 * (log(l_aa ^2 * l_ab ^2 + l_aa ^2 * l_bb ^2))
  ln_gam_bit <- -log(2 * sqrt(pi)) + lgamma((delta_star + numSq2 + 1) / 2 ) - lgamma((delta_star + numSq2 + 2) / 2)
  ln_h <- ln_top - ln_bottom + ln_gam_bit
  output <- ln_h
  
  return(output)
}

########## generate_graph_zo ##########

generate_graph_zo <- function(g_current, delta, Phi, N, S_N, warmup, iters_to_keep)
{
  p <- max(dim(Phi)[1], dim(Phi)[2])
  iter <- iters_to_keep
  if (iter <= warmup)
  {
    stop('not enough iterations specified')
  }
  sample <- iter - warmup
  counter_for_keeping <- 0
  graph_cumulative <- matrix(0, p, p)
  graph_pm <- matrix(0, p, p)
  graph_iter <- array(0, c(p, p, iters_to_keep))
  graph_tally <- array(0, c(p, p, 1))
  num_of_diff_graphs <- 1
  accept_count_warmup <- 0
  accept_count_sample <- 0
  accept_percent_warmup <- 0
  accept_percent_sample <- 0
  g_next <- g_current
  for (i in 1 : iter)
  {
    equal_indicator <- 0
    g_current <- g_next
    ss <- check_chordal(g_current)
    chordal <- ss[[1]]
    order <- ss[[2]]
    cliques <- chordal_to_ripcliques_zo(g_current, order)
    jtree <- ripcliques_to_jtree_zo(cliques)
    sepsize <- sepsize_zo(cliques, jtree)
    reach_graph <- reachability_graph(g_current)
    u <- runif(1)
    accept_prob <- 0
    output <- next_graph_candidate_zo(g_current, jtree, sepsize, cliques, reach_graph)
    g_proposal <- output[[1]]
    a <- output[[2]]
    b <- output[[3]]
    C <- output[[4]]
    CASE_add <- output[[5]]
    CASE_delete <- output[[6]]
    if (CASE_add == 1)
    {
      ln_likelihood_ratio <- ln_h_ratio_zo(C, a, b, delta, Phi) - ln_h_ratio_zo(C, a, b, delta + N, Phi + S_N)
      likelihood_ratio <- exp(ln_likelihood_ratio)
    }else{
      if (CASE_delete == 1)
      {
        ln_likelihood_ratio <- ln_h_ratio_zo(C, a, b, delta + N, Phi + S_N) - ln_h_ratio_zo(C, a, b, delta, Phi)
        likelihood_ratio <- exp(ln_likelihood_ratio)
      }
    }
    accept_prob <- min(1, likelihood_ratio)
    if (u <= accept_prob)
    {
      g_next <- g_proposal
      if (i <= warmup)
      {
        accept_count_warmup <- accept_count_warmup + 1
      }else{
        accept_count_sample <- accept_count_sample + 1
      }
    }else{
      g_next <- g_current
    }
  }
  
  accept_percent_sample <- accept_count_sample / iters_to_keep
  output <- list(g_next, accept_percent_sample)
  
  return(output)
}
