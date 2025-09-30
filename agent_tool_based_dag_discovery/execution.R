# executable file for the agent based LLM discovery

source("~/Desktop/master_thesis/code/MasterThesis/agent_tool_based_dag_discovery/tools_agent_LK.R")         # the tool implementations
source("~/Desktop/master_thesis/code/MasterThesis/agent_tool_based_dag_discovery/agent_llm_controller.R")    # the planner/controller

D_obs <- read_data(int = "none", path = "~/Desktop/master_thesis/code/MasterThesis/agent_tool_based_dag_discovery/data")

out <- run_llm_agent(D_obs, variables = colnames(D_obs), max_steps = 8, alpha = 0.05)

# Quick summary like before:
library(dplyr)
out$last_results %>%
  group_by(test) %>%
  summarise(n_cis_tested = unique(n_cis_tested),
            falsified = any(adj.p.value < 0.05),
            n_violations = sum(adj.p.value < 0.05))
