abstract type FiniteDifferencesMethod end

struct CentralFiniteDifferences <: FiniteDifferencesMethod
    order::Int
    derivative::Int
end

struct ForwardFiniteDifferences <: FiniteDifferencesMethod
    order::Int
    derivative::Int
end

struct BackwardFiniteDifferences <: FiniteDifferencesMethod
    order::Int
    derivative::Int
end
