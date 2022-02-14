from query import QueryNCBI

query = QueryNCBI(
    journal="Nature",
    year="2021"
    )
    
for res in query.results_iter():
    print(res)