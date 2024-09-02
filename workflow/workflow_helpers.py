def get_venues(config):
    return(list(config["venues"].keys()))

def get_novelty_null_iters(config):
    return list(range(0, config["novelty"]["iterations"]))

def get_novelty_years(config):
    return (
        list(
            range(
                config["novelty"]["start_year"],
                config["novelty"]["end_year"]
            )
        )
    )