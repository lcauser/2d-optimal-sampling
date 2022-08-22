#=
    Observers for the doob dynamics
=#

# Left components
function leftMeasure(model)
    return copy(model.leftComponent)
end

function leftObserver()
    return observer("left", "configuration", leftMeasure, 0, 0, 0)
end

# Escape rate
function escapeRateMeasure(model)
    return copy(sum(model.transitionRates))
end

function escapeRateObserver()
    return observer("ER", "configuration", escapeRateMeasure, 0, 0, 0)
end

# Original escape rate
function originalEscapeRateMeasure(model)
    return copy(sum(model.originalTransitionRates))
end

function originalEscapeRateObserver()
    return observer("OER", "configuration", originalEscapeRateMeasure, 0, 0, 0)
end
