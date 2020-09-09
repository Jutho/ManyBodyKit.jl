abstract type QuantumOperator{Q} end

abstract type LocalOperator{S} end <: QuantumOperator{S}
