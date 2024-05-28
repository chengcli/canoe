
option(FETCH_HITRAN "Fetch HITRAN data files" ON)
option(FETCH_JUPITER "Fetch Jupiter data files" ON)
option(FETCH_EXOGCM "Fetch ExoGCM data files" ON)

# Add a custom command that runs the script
add_custom_target(CheckDataIntegrity ALL
    COMMAND chmod +x ${CMAKE_CURRENT_SOURCE_DIR}/check_data_integrity.sh
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_data_integrity.sh ${CMAKE_CURRENT_SOURCE_DIR} ${FETCH_HITRAN} ${FETCH_JUPITER} ${FETCH_EXOGCM}
    COMMENT "Running check_data_integrity.sh"
    VERBATIM
)
