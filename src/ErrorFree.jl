@inline FastTwoSum(x, y) = (x+y, y-((x+y)-x))
@inline FastTwoDiff(x,y) = (x-y, (x-(x-y))-y)

@inline function TwoSum(x, y) 
  	s = x+y
  	aa = s-y
  	bb = s-aa
  	da = x-aa
  	db = y-bb
  	return s, da+db
end

@inline function TwoDiff(a, b)
	s = a-b
	bb = s-a
	aa = s-bb
	da = a-aa
	db = b+bb
	return s, da-db
end

@inline TwoProd(a, b) = (a*b, fma(a, b, -a*b))
