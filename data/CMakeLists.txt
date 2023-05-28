
# Add a custom command that runs the script
add_custom_target(CheckDataIntegrity ALL
    COMMAND chmod +x ${CMAKE_CURRENT_SOURCE_DIR}/check_data_integrity.sh
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_data_integrity.sh ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running check_data_integrity.sh"
    VERBATIM
)
