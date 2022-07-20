function strcoeff(args...)
    buf = IOBuffer()
    EnergyExpressions.showcoeff(buf, args...)
    String(take!(buf))
end

@testset "Pretty-print coefficients" begin
    @test strcoeff(1, false) == ""
    @test strcoeff(1, true) == "+ "
    @test strcoeff(-1, false) == "- "
    @test strcoeff(-1, true) == "- "

    @test strcoeff(2, false) == "2.0"
    @test strcoeff(2, true) == "+ 2.0"
    @test strcoeff(-2, false) == "- 2.0"
    @test strcoeff(-2, true) == "- 2.0"

    @test strcoeff(1, false, true) == "1"
    @test strcoeff(1, true, true) == "+ 1"
    @test strcoeff(-1, false, true) == "- 1"
    @test strcoeff(-1, true, true) == "- 1"

    @test strcoeff(im, false) == "im"
    @test strcoeff(im, true) == "+ im"
    @test strcoeff(-im, false) == "- im"
    @test strcoeff(-im, true) == "- im"

    @test strcoeff(2im, false) == "2.0im"
    @test strcoeff(2im, true) == "+ 2.0im"
    @test strcoeff(-2im, false) == "- 2.0im"
    @test strcoeff(-2im, true) == "- 2.0im"

    @test strcoeff((1+im), false) == "(1.0 + 1.0im)"
    @test strcoeff((1+im), true) == "+ (1.0 + 1.0im)"
    @test strcoeff(-(1+im), false) == "(-1.0 - 1.0im)"
    @test strcoeff(-(1+im), true) == "+ (-1.0 - 1.0im)"
end
