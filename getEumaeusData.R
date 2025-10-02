## 10/01/2025
## download some EUMAEUS example profile likelihoods for testing

library(dplyr)

# connection details
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))

# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

# get exposure information
sql = 'SELECT * from eumaeus.EXPOSURE'
exposures = DatabaseConnector::querySql(connection, sql)
names(exposures) = tolower(names(exposures))
print(exposures)
## can pick for case study: 21184 (H1N1), 21185 (Fluvirin), 21214 (Fluzone), 211983 (1/2 dose of Zoster), 211833 (1/2 dose of HPV)

# get list of negative control outcomes
sql = "SELECT * FROM eumaeus.negative_control_outcome"
NCs = DatabaseConnector::querySql(connection, sql)
names(NCs) = tolower(names(NCs))
NC_outcome_ids = NCs$outcome_id

# get some example profile likelihoods -- for NCs only!
sql <- "SELECT * FROM eumaeus.likelihood_profile
        WHERE method = 'HistoricalComparator'
              AND database_id = 'CCAE'
              AND analysis_id = 1
              AND exposure_id = 211983
              AND outcome_id IN (@outcomeIds);"
sql <- SqlRender::render(sql, outcomeIds = NC_outcome_ids)
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
exampleLPs <- DatabaseConnector::querySql(connection, sql)
names(exampleLPs) = tolower(names(exampleLPs))

glimpse(exampleLPs)
save(exampleLPs, file = "data/CCAE_HC_profiles.RData")

disconnect(connection)
