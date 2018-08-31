## Test history.stamp and a demonstration of provenance

test.history.stamp <- function(x) {
  attr(y,'history') <- history.stamp(x)
  return(y)
}

t2m <- t2m.DNMI()
y <- test.history.stamp(t2m)
print(provenance(y))